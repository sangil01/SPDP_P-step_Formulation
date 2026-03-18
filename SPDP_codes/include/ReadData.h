#ifndef SPDP_READ_DATA_H
#define SPDP_READ_DATA_H

#include <string>
#include <vector>

namespace spdp {

// 출발지, 도착지, 컨테이너 타입으로 구성된 하나의 요청 정보.
struct Request {
    int from_id = 0;
    int to_id = 0;
    int container_type = 0;
};

// `.dat` 파일에서 읽어온 SPDP 인스턴스 전체 데이터.
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

// 인스턴스 파일을 읽어 스칼라 값, 요청 목록, 행렬 정보를 채운다.
SPDPData read_spdp_data(const std::string& filename);

}  // namespace spdp

#endif
