import std.stdio, std.math, std.algorithm;
import std.array, std.container, std.random;


//To remove
//http://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm
//http://nghiaho.com/?page_id=355

mixin template MatrixT (T){
	//In case when squre matrix
	T[][] transpose(T [][] matrix)
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

	T[][] product(ref T[][] matrix, T number){
		return map!(x => map!(y => y * number)(x).array)(matrix).array;
	}

	T[][] product(T[][] matrix, T[][] values){
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

	T[] product(T[][] matrix, T[] vector){
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


	T[][] minor(T[][] matrix){
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
	T[][] adjoint(T[][] matrix){
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
	IMatrix jacobian(double input, double function(double, double[]) func){
		auto matrix = _data;
		for(int i = 0;i < matrix.length;++i){
			for(int j = 1;j <= matrix.length;++j)
				matrix[i][j] = Derivative(matrix[i], matrix[i], input, func, j);
		}
		return new Matrix(matrix);
	}
}


//http://en.wikipedia.org/wiki/Derivative
double Derivative(double []values1, double [] values2, double input, 
	double function(double, double[]) func, int size){

	double result1 = func(input, values1);
	double result2 = func(input, values2);
	return (result1 - result2)/size;
	//return (F(matrix) - F(matrix2Matrix))/matrix.length;
}


//TODO
double[][] JacobianResult(double[][]matrix, double inputvalue, double setvalue){
	mixin MatrixT!double;
	return product(transpose(matrix), matrix);
}



double GaussNewton(int iterations, double input[], double observed[], uint[string] data,
	double function(double, uint[string]) func)
	in{
		assert(input.length == observed.length);
	}
	body{
	int m = input.length;
	auto matr = new Matrix(m);
	double[][] Jacobian;
	double result = 0;
	double minerror = 10000000000;
	for(int i = 0;i < iterations;++i){
		auto param1 = new double[observed.length];
		double[] rdata = new double[m];
		double error = 0;
		for(int j = 0;j < m;++j){
			rdata[j] = observed[j] - func(input[j], data);
			error += (rdata[j] * rdata[j]);
			for(int k = 0;k < input.length;++k)
				auto newmatrix = matr.jacobian(input[k], func);
		}
		/*param1[i] = valueres;
		if(result < minerror){
			minerror = result;
		}
		auto value = matr * matr.T();
		auto value2 = value.inv();*/
		//need m-v product
		//auto ee = ((value2 * matr.T()) * value2) * param1;
	}
	return minerror;
}

void test_matrix(){
	//writeln(transpose([[1,2,3], [4,5,6], [7,8,9]]));
	auto data = [[1.0,2.0,3.0], [4.0,5.0,6.0], [7.0,8.0,9.0]];
	auto data2 = [[7.0,4.0, 5.0], 
				 [9.0,8.0,5.0], 
				 [2.0,1.0,1.0]];
	auto values = [[3.0,2.0,8.0], [4.0,5.0,1.0], [7.0,9.0,9.0]];

	auto m = new Matrix(data);
	auto m2 = m * data2;
	//writeln(m2.data);
}

double F(double inpvalue, double []values){
	double a = values[0];
	double b = values[1];
	double c = values[2];
	double d = values[3];
	return a * cos(b * inpvalue) + c * sin(d * inpvalue);
}

double targetFunc(double value, uint[string] data){
	uint a = data["A"];
	uint b = data["B"];
	uint c = data["C"];
	return a * cos(b * value) + sin(c * value);
}

double[] generateData(int count){
	double[]v = new double[count];
	return map!(x => uniform(-50.0, 50.0))(v).array;
}

double[] generateOutputdata(double function(double, uint[string]) func, uint[string] data, int count){
	double[]o = new double[count];
	for(int i = 0;i < count;++i){
		o[i] = func(i, data);
	}
	return o;
}

void test_gauss_newton(){
	uint[string] data;
	data["A"] = 5;
	data["B"] = 2;
	data["C"] = 7;
	double[] gendata = generateData(100);
	//auto otp = generateOutputdata((double x) => cast(double)cos(x),100);
	auto otp = generateOutputdata(&targetFunc, data, 100);
	auto result = GaussNewton(200, gendata, otp, data, &targetFunc);
}

void main()
{
	//writeln(F(0.4, [2.0,3.0,2.0,1.0]));
	test_gauss_newton();
}