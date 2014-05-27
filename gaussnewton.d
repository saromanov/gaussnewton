import std.stdio, std.math, std.algorithm;
import std.array, std.container;



mixin template Matrix (T){
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

	//Matrix inverse;
	//A^-1 = 1/|A|
	T [][] inverse(T [][] matrix){
		return matrix;
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

	// Get 2x2 matrix
	private double computeDet(T[][]matrix, T value)
	in{
		assert(matrix.length == 2 && matrix[0].length);
	  }
	body
	{
		return value * (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
	}
}

//Compute Jacobian Matrix
double[][] Jacobian(double[][]matrix, double input,double function(double, double[]) func)
{
	for(int i = 0;i < matrix.length;++i){
		for(int j = 1;j <= matrix.length;++j)
			matrix[i][j] = Derivative(matrix[i], matrix[i], input, func, j);
	}
	return matrix;
}

//http://en.wikipedia.org/wiki/Derivative
double Derivative(double []values1, double [] values2, double input, 
	double function(double, double[]) func, int size){

	double result1 = func(input, values1);
	double result2 = func(input, values2);
	return (result1 - result2)/size;
	//return (F(matrix) - F(matrix2))/matrix.length;
}


//TODO
double[][] JacobianResult(double[][]matrix, double inputvalue, double setvalue){
	mixin Matrix!double;
	return product(transpose(matrix), matrix);
}


double GaussNewton(int iterations, double input[], double observed[], double[][] matrix, 
	double function(double, double[]) func, double[] values)
	in{
		assert(input.length == observed.length);
	}
	body{

	double[][] Jacobian;
	int m = matrix.length;
	int n = matrix[0].length;
	double result = 0;
	for(int i = 0;i < iterations;++i){
		for(int j = 0;j < m;++j){
			//x_j + 1 = x_j - H * grad(f(x*))
			//double[][] newmatrix = Jacobian(matrix, input[i], func);
			double param1 = observed[j] - func(input[j], values);
			result += pow(param1,2);
		}
	}
	return result;
}

void test_matrix(){
	mixin Matrix!double;
	//writeln(transpose([[1,2,3], [4,5,6], [7,8,9]]));
	auto data = [[1.0,2.0,3.0], [4.0,5.0,6.0], [7.0,8.0,9.0]];
	auto data2 = [[7.0, 4.0, 5.0], [9.0,8.0,5.0], [2.0,1.0,1.0]];
	auto values = [[3.0,2.0,8.0], [4.0,5.0,1.0], [7.0,9.0,9.0]];
	auto prod = product(data, values);
	auto deter = determinant(data2);
	writeln(deter);
	//writeln(product(data, 5));
}

double F(double inpvalue, double []values){
	double a = values[0];
	double b = values[1];
	double c = values[2];
	double d = values[3];
	return a * cos(b * inpvalue) + c * sin(d * inpvalue);
}

void test_gauss_newton(){
	
}

void main()
{
	//writeln(F(0.4, [2.0,3.0,2.0,1.0]));
	test_matrix();
}