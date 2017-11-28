compile:
	g++ -o output main.cpp Individual.cpp IBEA.cpp WFG2.cpp WFG3.cpp WFG8.cpp -std=c++11
clean:
	rm output

exect:
	./output WFG2

# for i in `seq 0 4`; do ./hv ../results/WFG8_solution$i\_scale_0.005000.txt -r "3 5"; done
