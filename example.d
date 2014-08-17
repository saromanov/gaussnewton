import gaussnewton;
import std.stdio, std.math, std.algorithm;
import std.array, std.container, std.random, std.range;

double F(double inpvalue, double []values){
	double a = values[0];
	double b = values[1];
	double c = values[2];
	double d = values[3];
	return a * cos(b * inpvalue) + c * sin(d * inpvalue);
}

double targetFunc(double value, VARS data){
	float a = data["A"];
	float b = data["B"];
	float c = data["C"];
	return a * cos(b * value) + sin(c * value);
}

double[] generateData(int count){
	double[]v = new double[count];
	return map!(x => uniform(-50.0, 50.0))(v).array;
}


void test_gauss_newton(){
	auto data = [0.458,0.555,0.687];
	auto gn = new GaussNewton();
	gn.addVariable("A",5.0);
	gn.addVariable("B",3.0);
	gn.addVariable("C",8.0);
	gn.addFuncs(&targetFunc);
	gn.run(data,20);
}

void test_func_gauss_newton(){
	auto data = [0.333,0.5587,0.8225];
	auto gn = new GaussNewton();
	gn.addVariableF("A",0.85)
	  .addVariableF("B",0.99)
	  .addVariableF("C",0.22);
	gn.addFuncs(&targetFunc);
	writeln(gn.run(data, 20));
}

void test_gauss_newton_load(){
	auto gn = new GaussNewton();
	gn.fromFile("./data");
}


void main()
{
	test_func_gauss_newton();
	//test_gauss_newton_load();
}