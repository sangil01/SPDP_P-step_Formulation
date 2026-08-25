## 연구 개발

### Optimal travel cost, $c_{\min}$와 minimum vehicle number, $k_{\min}$을 동시에 달성하는 fractional solution을 자를 수 있는 기존 cut 반영 및 new cut 개발.

  - 새로운 problem-specific cut
  - Open-route VRP에서 사용되는 기존 cut 사용

### exact $k_{\min}$ 또는 더 엄밀한 $k_{\min}$ 계산

- C19와 D5의 경우 현재 구한  $k_{\min}$ 으로 incumbent infeasible을 600s 이내 판정
- C16, C18, C19, C20, D6은 현재 구한 $k_{\min}$으로 infeasibility 판정이 잘 되지 않음.
	1. 시간을 더 주고 infeasibility 판정이 된다면 $k_{\min} = k_{\min}+1$ 로 설정하여 풀면 됨. 
	2. CP로 feasibility check가 더 빠를 수 있음.
	3. DFF-based $k_{\min}$을 검토한다. 기존 duration bound는
	
	$$
	k_{\min}
	=
	\min_{y\in\mathcal Y}
	\left\lceil\frac{\sum_{e\in E}t_e y_e}{T}\right\rceil
	$$
	
	이고, DFF를 적용하면
	
	$$
	\min_{y\in\mathcal Y}
	\sum_{e\in E}f\!\left(\frac{t_e}{T}\right)y_e
	$$
	
	형태의 LP/IP를 고려할 수 있다. 이 접근은 유효성의 문제는 없지만, DFF를 고정된 request workload가 아니라 선택되는 multigraph edge duration에 적용한다는 한계가 있다. **DFF에 입력되는 edge 집합과 duration 분할 자체가 $y$에 따라 달라지므로, sub-IP가 DFF transformed value가 작은 edge들로 구성된 해를 선택할 수 있다. 따라서 route-level workload와 request 간 bin-packing incompatibility를 충분히 포착하지 못해 $k_{\min}$를 강화하지 못할 수 있다.**
	
- Duration IP를 풀 경우에는 parallel edge들 중 travel time이 가장 작은 edge를 제외하고는 추가 pruning을 한다. `empty state - routing - empty state` edge 또한 추가 pruning을 한다.
- Duration IP를 끝까지 풀 필요 없이 solve 도중 $\left\lceil\frac{UB}{T}\right\rceil=\left\lceil\frac{LB}{T}\right\rceil$ 가 성립하면 조기 종료한다.

### Algorithm 및 formulation 개발

- C8의 경우 duration IP based incumbent를 사용한 경우 optimal을 구했으나, 해당 OFV가 original travel-cost IP의 optimal solution을 warm start로 제공했을 때 구한 warm start OFV 보다 크다는 issue가 있어 이를 분석.
- Large-scale instance에는 PALNS 기반 two-stage heuristic 같은 방법을 검토하고 SOTA와 비교.
- Fixed-SPDP도 풀기 (Empty state를 type이 아니라 pickup location으로 가짐) → graph size가 커질 여지 있음. (state는 다양해지고, edge 수는 줄어들 수도)
- 동일 pick-up request (pickup location, treatment location, skip type 동일)/ 동일 delivery request (skip type, delivery location 동일) merge → disconnected cycles: root node & integer solution에서는 separation으로 cut 추가 (Forbes et al.)

### Generalization

- SPDP 뿐만 아니라, open route with fixed fleet cost VRP variants는 잘 풀릴 수 있음. exact가 잘 알려지지 않은 다른 VRP variant 찾아서 추가 실험.

## Forbes et al. 논문과의 비교

- Forbes et al.과 현재 연구의 핵심 기여 및 차이점을 명확히 구분한다.
	- Forbes et al.은 LP relaxation에서 차량 수를 먼저 최소화하고, 그 결과로 얻은 $k_{\min}$으로 차량 수를 고정한 뒤 travel cost를 최소화한다. 현재 연구는 계산한 $k_{\min}$을 valid inequality로 추가하고 $k_{\min}$를 고정하여 구한 incumbent solution을 original-cost IP의 warm start로 사용한다.
	- Forbes et al.의 $k_{\min}$이 우리 연구보다 좋지 않다.
	- Fragment-based graph와 우리의 multigraph가 유사하다.
		  - Forbes et al.: pickup 부터 delivery가 발생한 이후 발생한 pickup 직전까지의 restricted fragment를 하나의 arc로 표현한다. Fragment 내부에서는 pickup이 2번 가능하지만 모든 delivery보다 먼저 일어난다.
		  - 현재 연구: pickup/delivery request별 node를 만들고, 두 node 사이의 가능한 state pair 및 중간 emptying sequence를 parallel multigraph edge로 표현한다. 따라서 하나의 multigraph edge가 fragment보다 훨씬 수가 적다. (1)
- Forbes et al.은 pickup location, treatment location 및 skip type이 같은 request들을 multiplicity $q_r$를 갖는 하나의 request로 aggregate한다. 이 구조에서 발생하는 disconnected cycle을 separation한다. 그러나 우리 연구는 한 번만 방문을 강제하기에 subtour elimination을 flow-based constraint로 separation 없이 바로 넣을 수 있다. (2)
- (1), (2)에 의해 훨씬 빠르게 optimal을 구한다.
- 동일한 optimal instance를 기준으로 다음을 비교한다.
	  - root LP gap
	  - runtime
	  - graph 및 Gurobi model size: node 수, edge 수, 변수·제약 수
- LP gap은 유사하지만 model size는 작고 runtime이 빠른 것을 통해, multi graph가 route structure를 잘 반영하면서 model size도 작은 더 효율적인 graph라고 주장할 수 있다. 
