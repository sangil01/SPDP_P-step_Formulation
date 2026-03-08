#ifndef SPDP_READ_DATA_H
#define SPDP_READ_DATA_H

#include <string>
#include <vector>

namespace spdp {

struct Request {
    int from_id = 0;
    int to_id = 0;
    int container_type = 0;
};

struct SPDPData {
    double fixed_vehicle_cost = 0.0;
    double time_pickup = 0.0;
    double time_empty = 0.0;
    double time_delivery = 0.0;
    double time_limit = 0.0;
    int locations = 0;

    std::vector<Request> requests;
    std::vector<std::vector<double>> distance;
    std::vector<std::vector<double>> time;
};

SPDPData read_spdp_data(const std::string& filename);

}  // namespace spdp

#endif
