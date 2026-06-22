# Issue 1. P = 3일때 pricing이 P = 2일때 full enumeration 보다 좋지 않음

- Bucket graph 적용해볼 수 있는지 검토 (논문 읽고 세미나 준비까지)

# Issue 2. P = 2일 때 full enumeration을 한 게 flow based MILP보다 root node LP bound가 훨씬 작은 경우

1. 인스턴스 구조 분석
	- 특히, 다양한 구조의 인스턴스를 추가 생성하여 어떤 인스턴스가 문제인건지 확인
	
2. Gurobi 세팅 분석
3. ORtools CP로 풀어보기 (pattern-based model이 CP에 적합할 수 있음)
4. Full enumeration을 해놓고 Master에는 일부만 넣은 뒤 full column에 대해 RC 반복 update하는 BnP
	- RC 계산을 Matric form으로 할 수 있으면, GPU 병렬 컴퓨팅 가능
5. Heuristic solution을 p-step으로 쪼개서 초기 bound 주기 (warm start)
6. 큰 인스턴스도 풀어보기 (후순위)
