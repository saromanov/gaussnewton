module gaussnewton;

import std.stdio, std.math, std.algorithm;
import std.array, std.container, std.random, std.range;
import std.file;
import std.conv;

//http://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm


const float epsilon = 1e-4;
const float STEP = 1e-5;
alias float[string] VARS; //Variables
alias double function(double, float[string]) DFUNC;

mixin template MatrixT (T){
	alias T[][] MT; //Matrix type
	//In case when squre matrix
	MT transpose(T [][] matrix)
        in
        {
        	assert(matrix.length > 0 && matrix[0].length > 0 
        		&& matrix.length == matrix[0].length);
        }
        body {
		T [][] result = new T[][](matrix.length, matrix[0].length);
		for(int i = 0;i < matrix.length;++i){
			T tempvec[] = new T[matrix[i].length];
			for(int j = 0;j < matrix[i].length;++j){
				result[i][j] = matrix[j][i];
			}
		}
		return result;
	}

	MT product(ref T[][] matrix, T number){
		return map!(x => map!(y => y * number)(x).array)(matrix).array;
	}

	MT product(T[][] matrix, T[][] values){
		T [][] result = new T[][](matrix.length, matrix[0].length);
		for(int i = 0;i < matrix.length;++i){
			for(int j = 0;j < matrix[i].length;++j){
				double value = 0;
				for(int k = 0;k < values[j].length;++k)
					value += matrix[i][k] * values[k][j];
				result[i][j] = value;
			}
		}
		return result;
	}

	T[] product(MT matrix, T[] vector){
		T[] result = new T[](matrix.length);
		for(int i = 0;i < matrix.length;++i){
			T value = 0;
			for(int j = 0;j < matrix[0].length;++j){
				value += matrix[i][j] * vector[j];
			}
			result[i] = value;
		}
		return result;
	}

	//Matrix inverse;
	//A^-1 = 1/|A|
	//Now, only for 3x3 case
	T [][] inverse(T [][] matrix){
		T[][] adj = adjoint(matrix);
		return product(adj, (1/determinant(matrix)));
	}


	MT minor(T[][] matrix){
		T[][] result = new T[][](matrix.length, matrix.length);
		for(int p = 0;p < matrix.length;++p){
			for(int i = 0;i < matrix.length;++i){
				T [][]arr = new T[][](matrix.length-1, matrix.length-1);
				T [] temp;
				int a = 0;
				for(int j = 0;j < matrix.length;++j){
					for(int k = 0;k < matrix.length;++k){
						if(k == i)k+=1;
						if(j == p)j+=1;
						if(k == matrix.length)break;
						if(j == matrix.length)break;
						temp ~= matrix[j][k];
						if(temp.length == matrix.length-1){
							arr[a] = temp;
							temp = [];
							a+=1;
						}
					}
				}
				result[p][i] = computeMinor(arr);
			}

		}
		return transpose(result);
	}
    
    //http://en.wikipedia.org/wiki/Determinant
    //http://www.easycalculation.com/matrix/learn-matrix-determinant.php
	double determinant(T [][] matrix){
		double resultvalue = 0;
		for(int i = 0;i < matrix.length;++i){
			T value = matrix[0][i];
			T [][]arr = new T[][](matrix.length-1, matrix.length-1);
			T [] temp;
			int a = 0;
			for(int j = 1;j < matrix.length;++j){
				for(int k = 0;k < matrix.length;++k){
					if(k == i)k+=1;
					if(k == matrix.length)break;
					temp ~= matrix[j][k];
					if(temp.length == matrix.length-1){
						arr[a] = temp;
						temp = [];
						a+=1;
					}
				}
			}
		if(i == 1) resultvalue -= computeDet(arr, value);
		else
		resultvalue+=computeDet(arr, value);
		}
		return resultvalue;
	}

	//Find adjoint of a matrix where AB = I
	MT adjoint(T[][] matrix){
		T[][]min = minor(matrix);
		for(int i = 1;i <= min.length;++i){
			for(int j = 1;j <= min.length;++j){
				min[i-1][j-1] = pow(-1, i+j) * min[i-1][j-1];
			}
		}
		return min;
	}

	// Get 2x2 matrix
	double computeDet(T[][]matrix, T value)
	in{
		assert(matrix.length == 2 && matrix[0].length);
	  }
	body
	{
		return value * (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
	}

	double computeMinor(T[][] matrix)
	in{
		assert(matrix.length == 2 && matrix[0].length);
	  }
	body
	{
		return (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
	}
}

interface IMatrix
{
	IMatrix opMul(IMatrix);
	IMatrix opMul(double[][]);
	IMatrix opMul(double[]);
	IMatrix inv();
	IMatrix T();
	@property {
		double [][]data();
	}
}

class Matrix:IMatrix {
	private double[][]_data;
	this(double[][]data){
		_data = data;
	}

	this(double[] data){
		_data = new double[][](1, data.length);
		_data[0] = data;
	}

	this(int count){
		_data = new double[][](count, count);
	}

	@property {
		double [][]data(){ return _data;};
	}

	IMatrix opMul(IMatrix matr){
		return generalProduct(matr.data);
	}

	IMatrix opMul(double [][] matr){
		return generalProduct(matr);
	}

	IMatrix opMul(double[] vec){
		mixin MatrixT!double;
		return new Matrix(product(_data, vec));
	}

	private IMatrix generalProduct(double [][] matrix){
		mixin MatrixT!double;
		return new Matrix(product(_data, matrix));
	}

	IMatrix inv(){
		mixin MatrixT!double;
		return new Matrix(inverse(_data));
	}

	IMatrix T(){
		mixin MatrixT!double;
		return new Matrix(transpose(_data));
	}

	//Compute Jacobian of Matrix
	/*Matrix jacobian(double input, VARS data,
		double function(double, VARS data) func, int j){
		for(int i = 0;i < data.keys.length;++i){
			 data[j]  += 1;
			_data[i][j] = Derivative(input, func, data1, data2);
		}
		return new Matrix(_data);
	}*/
}


//http://en.wikipedia.org/wiki/Derivative
double Derivative(double[]input, 
	double function(double, VARS) func, VARS data1, VARS data2){
	double result1 = func(input[0], data1);
	double result2 = func(input[0], data2);
	return (result1 - result2)/epsilon;
}



void changeParameters(VARS params){

}


//TODO
double[][] JacobianResult(double[][]matrix, double inputvalue, double setvalue){
	mixin MatrixT!double;
	return product(transpose(matrix), matrix);
}



double[] GaussNewtonImpl(int iterations, double input[], double observed[], float[string] data,
	double function(double, float[string]) func, float step)
	in{
		assert(input.length == observed.length);
	}
	body{
	int m = input.length;
	auto matr = new Matrix(m);
	double[][] jacobi = new double[][](m, data.keys.length);
	double result = 0;
	double minerror = 10000000000;
	double[][]bvalues = new double[][](iterations+1, iterations+1);
	bvalues[0] = [0.85,0.50,-0.30];
	for(int i = 0;i < iterations;++i){
		auto param1 = new double[observed.length];
		double[] rdata = new double[m];
		double error = 0;
		auto vars1 = data.dup;
		for(int j = 0;j < m;++j){
			rdata[j] = observed[j] - func(input[j], data);
			error += (rdata[j] * rdata[j]);
			int c = 0;
			foreach(key;data.keys){
			 	vars1[key]  += STEP;
				jacobi[j][c] = Derivative(input, func, vars1, data);
				c += 1;
			}
		}

		matr = new Matrix(jacobi);
		if(error < minerror){
			minerror = error;
		}
		auto value = matr * matr.T();
		auto value2 = value.inv();

		//Parameters for update
		//TODO fix m-v product
		bvalues[i+1] = add(bvalues[i], ((((value2 * matr.T()) * value2) * rdata).data[0]));
	}
	return bvalues[0];
}

double[] add(double[]data1, double[]data2){
	double[] res = new double[](data1.length);
	for(int i = 0;i < data1.length;++i)
		res[i] = data1[i] + data2[i];
	return res;
}


double[] generateOutputdata(double function(double, VARS) func, VARS data, int count){
	double[]o = new double[count];
	for(int i = 0;i < count;++i){
		o[i] = func(i, data);
	}
	return o;
}


interface IGaussNewton{
	void addVariable(string name, float value);
	void addFuncs(DFUNC func);
	IGaussNewton addVariableF(string name, float value);
	string[] showVariables();
	double[] run(double[]data, int iters);
}

class GaussNewton:IGaussNewton {
	private VARS variables;
	DFUNC _func;
	this(){
	}

	//Append variable
	this(VARS variables2){
		variables = variables2;
	}

	//Append variable for compute target function.
	//Of course, depends of numbers on var-s on this function
	void addVariable(string name, float value){
		variables[name] = value;
	}

	IGaussNewton addVariableF(string name, float value){
		variables[name] = value;
		return new GaussNewton(variables);
	}

	//Output format is variable-value
	void fromFile(string filename){
		auto lines = cast(char[]) read(filename);
		foreach(ref str; lines.split("\n")){
			auto value = str.split();
			this.addVariable(to!string(value[0]), to!float(value[1]));
		}

	}

	//Add target function
	void addFuncs(DFUNC func){
		_func = func;
	}

	//Show variables which was add on addVariables
	string[] showVariables(){
		return variables.keys;
	}


	//"heart" of this class. Run G-N algorithm
	double[] run(double[]data, int iters){
		if(_func != null && variables.length > 0){
			auto gens = generateOutputdata(_func, variables, 3);
			return GaussNewtonImpl(iters, data, gens, variables, _func, epsilon);
		}
		return new double[](1);
	}
}


