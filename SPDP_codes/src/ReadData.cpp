#include "ReadData.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace spdp {
namespace {

namespace fs = std::filesystem;

enum class ParseSection {
    Header,
    Requests,
    Distance,
    Time,
};

std::string trim(const std::string& value) {
    std::size_t begin = 0;
    while (begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin])) != 0) {
        ++begin;
    }

    std::size_t end = value.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(value[end - 1])) != 0) {
        --end;
    }

    return value.substr(begin, end - begin);
}

std::vector<std::string> split_tokens(const std::string& line) {
    std::istringstream input(line);
    std::vector<std::string> tokens;
    std::string token;
    while (input >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

std::string to_lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

std::string to_upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::toupper(ch));
    });
    return value;
}

std::string canonicalize_key(const std::string& value) {
    std::string normalized;
    normalized.reserve(value.size());

    for (unsigned char ch : value) {
        if (std::isalnum(ch) != 0) {
            normalized.push_back(static_cast<char>(std::tolower(ch)));
        }
    }

    return normalized;
}

double parse_number(const std::string& token) {
    try {
        return std::stod(token);
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid numeric token: " + token);
    }
}

int parse_int_token(const std::string& token) {
    return static_cast<int>(parse_number(token));
}

fs::path find_project_root() {
    std::vector<fs::path> starts;

    try {
        starts.push_back(fs::weakly_canonical(fs::path(__FILE__)).parent_path());
    } catch (const std::exception&) {
        // Fall through and use cwd only.
    }

    starts.push_back(fs::current_path());

    for (const fs::path& start : starts) {
        fs::path current = start;
        while (!current.empty()) {
            const fs::path candidate_data_dir = current / "SPDP_data";
            if (fs::exists(candidate_data_dir) && fs::is_directory(candidate_data_dir)) {
                return current;
            }

            const fs::path parent = current.parent_path();
            if (parent == current) {
                break;
            }
            current = parent;
        }
    }

    return fs::path{};
}

fs::path resolve_data_path(const std::string& filename) {
    const fs::path input_path(filename);
    if (fs::exists(input_path) && fs::is_regular_file(input_path)) {
        return fs::absolute(input_path);
    }

    const fs::path project_root = find_project_root();

    std::vector<fs::path> candidates;
    if (!project_root.empty()) {
        candidates.push_back(project_root / "SPDP_data" / filename);
    }
    candidates.push_back(fs::current_path() / "SPDP_data" / filename);
    candidates.push_back(fs::current_path() / filename);

    for (const fs::path& candidate : candidates) {
        if (fs::exists(candidate) && fs::is_regular_file(candidate)) {
            return fs::absolute(candidate);
        }
    }

    throw std::runtime_error("Could not find data file: " + filename);
}

std::pair<std::vector<std::vector<double>>, std::size_t> parse_matrix(
    const std::vector<std::string>& lines,
    std::size_t start_idx,
    int locations,
    const std::string& section_name
) {
    std::vector<std::vector<double>> rows;
    std::size_t idx = start_idx;

    while (idx < lines.size() && static_cast<int>(rows.size()) < locations) {
        const std::string line = trim(lines[idx]);
        if (line.empty()) {
            ++idx;
            continue;
        }

        const std::vector<std::string> tokens = split_tokens(line);
        std::vector<double> row;
        row.reserve(tokens.size());

        for (const std::string& token : tokens) {
            try {
                row.push_back(std::stod(token));
            } catch (const std::exception&) {
                throw std::runtime_error(
                    "Non-numeric token found while reading " + section_name +
                    " matrix at line " + std::to_string(idx + 1) + ": " + line
                );
            }
        }

        if (static_cast<int>(row.size()) != locations) {
            throw std::runtime_error(
                section_name + " matrix row length mismatch at line " + std::to_string(idx + 1) +
                ": expected " + std::to_string(locations) +
                ", got " + std::to_string(row.size())
            );
        }

        rows.push_back(std::move(row));
        ++idx;
    }

    if (static_cast<int>(rows.size()) != locations) {
        throw std::runtime_error(
            section_name + " matrix row count mismatch: expected " + std::to_string(locations) +
            ", got " + std::to_string(rows.size())
        );
    }

    return {rows, idx};
}

std::vector<std::vector<double>> augment_metric_with_virtual_zero(
    const std::vector<std::vector<double>>& metric
) {
    const std::size_t size = metric.size();
    std::vector<std::vector<double>> augmented(size + 1, std::vector<double>(size + 1, 0.0));

    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = 0; j < size; ++j) {
            augmented[i][j] = metric[i][j];
        }
    }

    return augmented;
}

}  // namespace

SPDPData read_spdp_data(const std::string& filename) {
    const fs::path data_path = resolve_data_path(filename);

    std::ifstream input(data_path);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to open data file: " + data_path.string());
    }

    std::vector<std::string> raw_lines;
    std::string line;
    while (std::getline(input, line)) {
        raw_lines.push_back(line);
    }

    std::unordered_map<std::string, std::optional<double>> params = {
        {"fixedvehiclecost", std::nullopt},
        {"timepickup", std::nullopt},
        {"timeempty", std::nullopt},
        {"timedelivery", std::nullopt},
        {"timelimit", std::nullopt},
        {"locations", std::nullopt},
    };

    std::vector<Request> requests;
    std::optional<int> request_count;
    std::vector<std::vector<double>> distance;
    std::vector<std::vector<double>> time;

    std::size_t idx = 0;
    ParseSection section = ParseSection::Header;

    while (idx < raw_lines.size()) {
        const std::string stripped_line = trim(raw_lines[idx]);
        if (stripped_line.empty()) {
            ++idx;
            continue;
        }

        const std::vector<std::string> tokens = split_tokens(stripped_line);
        if (tokens.empty()) {
            ++idx;
            continue;
        }

        const std::string key = tokens[0];
        const std::string key_lower = to_lower(key);
        const std::string canonical_key = canonicalize_key(key);

        if (key_lower == "requests") {
            section = ParseSection::Requests;
            if (tokens.size() > 1) {
                request_count = parse_int_token(tokens[1]);
            }
            ++idx;
            continue;
        }

        if (canonical_key == "distance" || canonical_key == "distances") {
            if (!params["locations"].has_value()) {
                throw std::runtime_error("LOCATIONS must be parsed before Distance matrix.");
            }

            const int locations = static_cast<int>(params["locations"].value());
            auto parsed = parse_matrix(raw_lines, idx + 1, locations, "Distance");
            distance = std::move(parsed.first);
            idx = parsed.second;
            section = ParseSection::Distance;
            continue;
        }

        if (canonical_key == "time" || canonical_key == "times") {
            if (!params["locations"].has_value()) {
                throw std::runtime_error("LOCATIONS must be parsed before Time matrix.");
            }

            const int locations = static_cast<int>(params["locations"].value());
            auto parsed = parse_matrix(raw_lines, idx + 1, locations, "Time");
            time = std::move(parsed.first);
            idx = parsed.second;
            section = ParseSection::Time;
            continue;
        }

        if (section == ParseSection::Header) {
            const auto found = params.find(canonical_key);
            if (found != params.end() && tokens.size() > 1) {
                found->second = parse_number(tokens[1]);
            }
        } else if (section == ParseSection::Requests) {
            if (to_upper(key) == "ID") {
                ++idx;
                continue;
            }

            if (tokens.size() < 6) {
                throw std::runtime_error(
                    "Invalid request row format at line " + std::to_string(idx + 1) +
                    ": " + stripped_line
                );
            }

            Request request;
            request.from_id = parse_int_token(tokens[2]);
            request.container_type = parse_int_token(tokens[3]);
            request.to_id = parse_int_token(tokens[5]);
            requests.push_back(request);
        }

        ++idx;
    }

    std::vector<std::string> missing_params;
    missing_params.reserve(params.size());
    for (const auto& entry : params) {
        if (!entry.second.has_value()) {
            missing_params.push_back(entry.first);
        }
    }

    if (!missing_params.empty()) {
        std::sort(missing_params.begin(), missing_params.end());
        std::string message = "Missing required global parameters:";
        for (const std::string& key : missing_params) {
            message += " " + key;
        }
        throw std::runtime_error(message);
    }

    if (distance.empty()) {
        throw std::runtime_error("Distance matrix was not parsed.");
    }

    if (time.empty()) {
        throw std::runtime_error("Time matrix was not parsed.");
    }

    if (request_count.has_value() && request_count.value() != static_cast<int>(requests.size())) {
        throw std::runtime_error(
            "REQUESTS count mismatch: header=" + std::to_string(request_count.value()) +
            ", parsed=" + std::to_string(requests.size())
        );
    }

    SPDPData data;
    data.fixed_vehicle_cost = params["fixedvehiclecost"].value();
    data.time_pickup = params["timepickup"].value();
    data.time_empty = params["timeempty"].value();
    data.time_delivery = params["timedelivery"].value();
    data.time_limit = params["timelimit"].value();
    data.locations = static_cast<int>(params["locations"].value());
    data.requests = std::move(requests);
    data.distance = augment_metric_with_virtual_zero(distance);
    data.time = augment_metric_with_virtual_zero(time);

    return data;
}

}  // namespace spdp
