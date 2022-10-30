#include <iostream>
#include <iomanip>
using namespace std ;
//在屏幕上打印矩阵
void PrintMat( double *Mat , int m , int n )
{
	int i , j ;
	for( i = 0 ; i < m ;i ++ )
	{
		for( j = 0; j < n ; j ++ )
			cout<<setw(15)<<Mat[ i * n + j ];
		cout<<endl;
	}
}

//矩阵置零
void MatBeZero( double *Mat , int m , int n )
{
	int i , j ;
	for( i = 0 ; i < m ; i ++ )
	for( j = 0 ; j < n ; j ++ )
		Mat[ i * n + j ] = 0.0 ;
}

//Mt(nxm) = transpose( M(mxn) )
void Matrix_Transpose( double *Mt , double *M , int m , int n )
{
	int i , j ;
	for( i = 0 ; i < n ; i ++ )
	{
		for( j = 0 ; j < m ; j ++ )
			Mt[ i * m + j ] = M[ j * n + i ] ;//Mt[i][j]=M[j][i]
	}
}
//A = B + C
void MatAdd( double *A , double *B , double *C , int m , int n )
{
	int i , j ;
	for( i = 0 ; i < m ; i ++ )
	for( j = 0 ; j < n ; j ++ )
		A[ i * n + j ] = B[ i * n + j ] + C[ i * n + j ] ;
}
//[ A ] = [ B ] * c
void MatMul( double *A , double *B , double c , int m , int n )
{
	int i , j ;
	for( i = 0 ; i < m ; i ++ )
	for( j = 0 ; j < n ; j ++ )
		A[ i * n + j ] = B[ i * n + j ] * c ;
}

//A(mxn)=B(mxp)*C(pxn)
void MatMul(double *A , double *B , double *C ,	
			int m , int p , int n)
{ 
	for(int i=0;i<m;i++)
	{
		for (int j=0;j<n;j++)
		{   
			A[i*n+j]=0;
			for(int k=0;k<p;k++)
			{
				A[i*n+j]+=B[i*p+k]*C[k*n+j];
			}
		}
	}
}

//基本高斯消元法
//A * x = y , A矩阵是m行*m列的矩阵
void Gauss( double *A , double *x , double *y , int m )
{
	int i , j , k ;
	double t ;//临时变量

	for( k = 0 ; k <= m - 2 ; k ++ )
	{
		//开始第k步消元，处理从第k行k列开始的子矩阵
		for( i = k + 1 ; i < m ; i ++ )
		{
			t = - A[ i * m + k ] / A[ k * m + k ] ;
			for( j = k ; j < m ; j ++ )
			{
				A[ i * m + j ] = A[ i * m + j ] + t * A[ k * m + j ] ;	
			}
			y[ i ] = y[ i ] + t * y[ k ] ;
		}
	}
	//回代求解
	x[ m - 1 ] = y[ m - 1 ] / A[ ( m - 1 ) * m + ( m - 1 ) ] ;
	for( i = m - 2 ; i >= 0 ; i -- )
	{
		t = y[ i ] ;
		for( j = i + 1 ; j < m ; j ++ )
			t = t - A[ i * m + j ] * x[ j ] ;
		x[ i ] = t / A[ i * m + i ] ;
	}
}

int main()
{
	double E , I ,  A , P ;//E:弹模 , I:惯性矩 , A:梁横截面积, P:集中载荷 
	double L0 , L1 ,L2;//单元0,1,2的长度
	double *ke0 , *ke1,*ke2;//单刚阵
	double alpha,beta ;//坐标转换角
	double T[ 36 ] , Tt[ 36 ] ;//坐标转换矩阵
	double *K ;//总刚阵
	double Tem[ 36 ] ;//临时矩阵
	int LOCEL0[ 6 ] , LOCEL1[ 6 ] ,LOCEL2[6];//单元0,1,2的门票
	int i , j , ii , jj ;
	double *F ,*F0,*F1,*F2, *V,*V0,*V1,*V2 ;//载荷向量, 位移向量

	E = 2.1e6 ;//弹模 kgf/cm^2
	I = 1.0 / 12;//惯性矩cm^4
	A = 1 ;//梁横截面积cm^2
	L0 = 200 ;//单元0长度200cm
	L1 = 200 * sqrt(2.0) ;//单元1长度
	L2 = 200 * sqrt(2.0);//单元2长度
	P = -10 ;//P=10 kgf
	//申请存储空间
	ke0 = new double[ 6 * 6 ] ;
	ke1 = new double[ 6 * 6 ] ;
	ke2 = new double[6 * 6];
	K = new double[ 12 * 12 ] ;
	F = new double[ 12 ] ;
	F0 = new double[6];
	F1 = new double[6];
	F2 = new double[6];
	V = new double[ 12 ] ;
	V0 = new double[6];
	V1 = new double[6];
	V2 = new double[6];

	//计算单刚
	MatBeZero( ke0 , 6 , 6 ) ;//单刚置零
	ke0[ 0 * 6 + 0 ] = A ;	ke0[ 0 * 6 + 3 ] = -A ;
	ke0[ 1 * 6 + 1 ] = 12 * I / L0 / L0 ;		ke0[ 1 * 6 + 2 ] = 6 * I / L0 ;		ke0[ 1 * 6 + 4 ] = - 12 * I / L0 / L0 ;		ke0[ 1 * 6 + 5 ] = 6 * I / L0 ;
	ke0[ 2 * 6 + 1 ] = 6 * I / L0 ;				ke0[ 2 * 6 + 2 ] = 4 * I ;				ke0[ 2 * 6 + 4 ] = - 6 * I / L0 ;					ke0[ 2 * 6 + 5 ] = 2 * I;
	ke0[ 3 * 6 + 0 ] = - A ;	ke0[ 3 * 6 + 3 ] = A ;
	ke0[ 4 * 6 + 1 ] = - 12 * I / L0 / L0 ;	ke0[ 4 * 6 + 2 ] = - 6 * I / L0 ;		ke0[ 4 * 6 + 4 ] = 12 * I / L0 / L0 ;			ke0[ 4 * 6 + 5 ] = - 6 * I / L0 ;
	ke0[ 5 * 6 + 1 ] = 6 * I / L0 ;				ke0[ 5 * 6 + 2 ] = 2 * I  ;				ke0[ 5 * 6 + 4 ] = - 6 * I  / L0 ;				ke0[ 5 * 6 + 5 ] = 4 * I  ;
	for( i = 0 ; i < 6 ; i ++ )
	{
		for( j = 0 ; j < 6 ; j ++ )
		{
			ke0[ i * 6 + j ] *= E / L0 ;
		}
	}
	MatBeZero( ke1 , 6 , 6 ) ;//单刚置零
	ke1[ 0 * 6 + 0 ] = A ;	ke1[ 0 * 6 + 3 ] = -A ;
	ke1[ 1 * 6 + 1 ] = 12 * I / L1 / L1 ;		ke1[ 1 * 6 + 2 ] = 6 * I / L1 ;		ke1[ 1 * 6 + 4 ] = - 12 * I / L1 / L1 ;		ke1[ 1 * 6 + 5 ] = 6 * I / L1 ;
	ke1[ 2 * 6 + 1 ] = 6 * I / L1 ;				ke1[ 2 * 6 + 2 ] = 4 * I ;				ke1[ 2 * 6 + 4 ] = - 6 * I / L1 ;					ke1[ 2 * 6 + 5 ] = 2 * I;
	ke1[ 3 * 6 + 0 ] = - A ;	ke1[ 3 * 6 + 3 ] = A ;
	ke1[ 4 * 6 + 1 ] = - 12 * I / L1 / L1 ;	ke1[ 4 * 6 + 2 ] = - 6 * I / L1 ;		ke1[ 4 * 6 + 4 ] = 12 * I / L1 / L1 ;			ke1[ 4 * 6 + 5 ] = - 6 * I / L1 ;
	ke1[ 5 * 6 + 1 ] = 6 * I / L1 ;				ke1[ 5 * 6 + 2 ] = 2 * I  ;				ke1[ 5 * 6 + 4 ] = - 6 * I  / L1 ;				ke1[ 5 * 6 + 5 ] = 4 * I  ;
	for( i = 0 ; i < 6 ; i ++ )
	{
		for( j = 0 ; j < 6 ; j ++ )
		{
			ke1[ i * 6 + j ] *= E / L1 ;
		}
	}
	MatBeZero(ke2, 6, 6);//单刚置零
	ke2[0 * 6 + 0] = A;	ke2[0 * 6 + 3] = -A;
	ke2[1 * 6 + 1] = 12 * I / L2/ L2;		ke2[1 * 6 + 2] = 6 * I / L2;		ke2[1 * 6 + 4] = -12 * I / L2 / L2;		ke2[1 * 6 + 5] = 6 * I / L2;
	ke2[2 * 6 + 1] = 6 * I / L2;				ke2[2 * 6 + 2] = 4 * I;				ke2[2 * 6 + 4] = -6 * I / L2;					ke2[2 * 6 + 5] = 2 * I;
	ke2[3 * 6 + 0] = -A;	ke2[3 * 6 + 3] = A;
	ke2[4 * 6 + 1] = -12 * I / L2 / L2;	ke2[4 * 6 + 2] = -6 * I / L2;		ke2[4 * 6 + 4] = 12 * I / L2 / L2;			ke2[4 * 6 + 5] = -6 * I / L2;
	ke2[5 * 6 + 1] = 6 * I / L2;				ke2[5 * 6 + 2] = 2 * I;				ke2[5 * 6 + 4] = -6 * I / L2;				ke2[5 * 6 + 5] = 4 * I;
	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < 6; j++)
		{
			ke2[i * 6 + j] *= E / L2;
		}
	}
	//
	cout<<"ke1="<<endl ;
	PrintMat( ke1 , 6 , 6 ) ;
	//1号单刚坐标转换
	alpha = 5.0 / 4 * acos(-1.0) ;//坐标转换角=5/4*PI
	MatBeZero( T , 6 , 6 ) ;
	T[ 0 * 6 + 0 ] = cos(alpha ) ;	T[ 0 * 6 + 1 ] = - sin( alpha ) ;
	T[ 1 * 6 + 0 ] = sin(alpha ) ;	T[ 1 * 6 + 1 ] =   cos( alpha ) ;
	T[ 2 * 6 + 2 ] = 1.0 ;
	T[ 3 * 6 + 3 ] = cos(alpha ) ;	T[ 3 * 6 + 4 ] = - sin( alpha ) ;
	T[ 4 * 6 + 3 ] = sin(alpha ) ;	T[ 4 * 6 + 4 ] =   cos( alpha ) ;
	T[ 5 * 6 + 5 ] = 1.0 ;
	Matrix_Transpose( Tt , T , 6 , 6 ) ;
	//
	cout<<"T="<<endl ;
	PrintMat( T , 6 , 6 ) ;
	//Tem = T * ke1 ;
	MatMul( Tem , T , ke1 , 6 , 6 , 6 ) ;
	//ke1 = Tem * Tt ;
	MatMul( ke1 , Tem , Tt , 6 , 6 , 6 ) ;
	//
	cout<<"ke1="<<endl ;
	PrintMat( ke1 , 6 , 6 ) ;
	//2号单刚坐标转换
	beta = 3.0 / 4 * acos(-1.0);//坐标转换角=3/4*PI
	MatBeZero(T, 6, 6);
	T[0 * 6 + 0] = cos(beta);	T[0 * 6 + 1] = -sin(beta);
	T[1 * 6 + 0] = sin(beta);	T[1 * 6 + 1] = cos(beta);
	T[2 * 6 + 2] = 1.0;
	T[3 * 6 + 3] = cos(beta);	T[3 * 6 + 4] = -sin(beta);
	T[4 * 6 + 3] = sin(beta);	T[4 * 6 + 4] = cos(beta);
	T[5 * 6 + 5] = 1.0;
	Matrix_Transpose(Tt, T, 6, 6);
	//
	cout << "T=" << endl;
	PrintMat(T, 6, 6);
	//Tem = T * ke1 ;
	MatMul(Tem, T, ke2, 6, 6, 6);
	//ke1 = Tem * Tt ;
	MatMul(ke2, Tem, Tt, 6, 6, 6);
	//
	cout << "ke2=" << endl;
	PrintMat(ke1, 6, 6);
	//组集总刚
	MatBeZero( K , 12 , 12 ) ;//总刚阵置零
	//0号单刚组集
	//单元0的门票
	LOCEL0[ 0 ] = 0 ;	LOCEL0[ 1 ] = 1 ;	LOCEL0[ 2 ] = 2 ;	LOCEL0[ 3 ] = 3 ;	LOCEL0[ 4 ] = 4 ;	LOCEL0[ 5 ] = 5 ;
	for( i = 0 ; i < 6 ; i ++ )
	{
		ii = LOCEL0[ i ] ;
		for( j = 0 ; j < 6 ; j ++ )
		{
			jj = LOCEL0[ j ] ;
			K[ ii * 12 + jj ] += ke0[ i * 6 + j ] ;
		}
	}
	//1号单刚组集
	//单元1的门票
	LOCEL1[ 0 ] = 3 ;	LOCEL1[ 1 ] = 4 ;	LOCEL1[ 2 ] = 5 ;	LOCEL1[ 3 ] = 6 ;	LOCEL1[ 4 ] = 7 ;	LOCEL1[ 5 ] = 8 ;
	for( i = 0 ; i < 6 ; i ++ )
	{
		ii = LOCEL1[ i ] ;
		for( j = 0 ; j < 6 ; j ++ )
		{
			jj = LOCEL1[ j ] ;
			K[ ii * 12+ jj ] += ke1[ i * 6 + j ] ;
		}
	}
	//2号单刚组集
	//单元2的门票
	LOCEL2[0] = 3;	LOCEL2[1] = 4;	LOCEL2[2] = 5;	LOCEL2[3] = 9;	LOCEL2[4] = 10;	LOCEL2[5] = 11;
	for (i = 0; i < 6; i++)
	{
		ii = LOCEL2[i];
		for (j = 0; j < 6; j++)
		{
			jj = LOCEL2[j];
			K[ii * 12 + jj] += ke2[i * 6 + j];
		}
	}

	//打印总刚
	cout<<"总刚:"<<endl ;
	PrintMat( K , 12 , 12 ) ;
	//载荷向量
	MatBeZero( F , 12 , 1 ) ;
	F[ 4 ] = P ;
	//边界条件处理
	//0号节点刚固(对角元置大数法)
	K[ 0 * 12 + 0 ] = 1e20 ;
	K[ 1 * 12 + 1 ] = 1e20 ;
	K[ 2 * 12 + 2 ] = 1e20 ;
	//2号节点刚固(对角元置大数法)
	K[ 6 * 12 + 6 ] = 1e20 ;
	K[ 7 * 12 + 7 ] = 1e20 ;
	K[ 8 * 12 + 8 ] = 1e20 ;
	//3号节点刚固(对角元置大数法)
	K[9 * 12 + 9] = 1e20;
	K[10 * 12 + 10] = 1e20;
	K[11 * 12 + 11] = 1e20;
	//求解
	Gauss( K , V , F ,  12) ;
	//打印节点位移向量
	cout<<"节点位移向量:"<<endl;
	PrintMat( V , 12 , 1 ) ;
	//求解0号单元杆端力
	MatBeZero(F0, 6, 1);
	MatBeZero(V0, 6, 1);
	V0[0 * 1 + 0] = V[0 * 12 + 0]; V0[1 * 1 + 0] = V[1 * 1 + 0]; V0[2 * 1 + 0] = V[2 * 1 + 0];
	V0[3 * 1 + 0] = V[3 * 1 + 0]; V0[4 * 1 + 0] = V[4 * 1 + 0]; V0[5 * 1 + 0] = V[5 * 1 + 0];
	//F0=Ke0*V0
	MatMul(F0, ke0, V0, 6, 6, 1);
	cout << "0号单元杆端力:" << endl;
	PrintMat(F0, 6, 1);
	//求解1号单元杆端力
	MatBeZero(F1, 6, 1);
	MatBeZero(V1, 6, 1);
	V1[0 * 1 + 0] = V[3* 12 + 0]; V1[1 * 1 + 0] = V[4 * 1 + 0]; V1[2 * 1 + 0] = V[5 * 1 + 0];
	V1[3 * 1 + 0] = V[6 * 1 + 0]; V1[4 * 1 + 0] = V[7 * 1 + 0]; V1[5 * 1 + 0] = V[8 * 1 + 0];
	//F1=Ke1*V1
	MatMul(F1, ke1, V1, 6, 6, 1);
	cout << "1号单元杆端力:" << endl;
	PrintMat(F1, 6, 1);
	//求解2号单元杆端力
	MatBeZero(F2, 6, 1);
	MatBeZero(V2, 6, 1);
	V2[0 * 1 + 0] = V[3 * 12 + 0]; V2[1 * 1 + 0] = V[4 * 1 + 0]; V2[2 * 1 + 0] = V[5 * 1 + 0];
	V2[3 * 1 + 0] = V[9 * 1 + 0]; V2[4 * 1 + 0] = V[10* 1 + 0]; V2[5 * 1 + 0] = V[11 * 1 + 0];
	//F2=Ke2*V2
	MatMul(F2, ke2, V2, 6, 6, 1);
	cout << "2号单元杆端力:" << endl;
	PrintMat(F2, 6, 1);
	
	//释放占用空间
	delete []ke0;
	delete []ke1;
	delete []ke2;
	delete []K ;
	delete []F ;
	delete []V ;

	system("Pause");
	return 0;
}