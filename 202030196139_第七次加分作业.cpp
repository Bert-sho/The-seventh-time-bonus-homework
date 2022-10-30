#include <iostream>
#include <iomanip>
using namespace std ;
//����Ļ�ϴ�ӡ����
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

//��������
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

//������˹��Ԫ��
//A * x = y , A������m��*m�еľ���
void Gauss( double *A , double *x , double *y , int m )
{
	int i , j , k ;
	double t ;//��ʱ����

	for( k = 0 ; k <= m - 2 ; k ++ )
	{
		//��ʼ��k����Ԫ������ӵ�k��k�п�ʼ���Ӿ���
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
	//�ش����
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
	double E , I ,  A , P ;//E:��ģ , I:���Ծ� , A:��������, P:�����غ� 
	double L0 , L1 ,L2;//��Ԫ0,1,2�ĳ���
	double *ke0 , *ke1,*ke2;//������
	double alpha,beta ;//����ת����
	double T[ 36 ] , Tt[ 36 ] ;//����ת������
	double *K ;//�ܸ���
	double Tem[ 36 ] ;//��ʱ����
	int LOCEL0[ 6 ] , LOCEL1[ 6 ] ,LOCEL2[6];//��Ԫ0,1,2����Ʊ
	int i , j , ii , jj ;
	double *F ,*F0,*F1,*F2, *V,*V0,*V1,*V2 ;//�غ�����, λ������

	E = 2.1e6 ;//��ģ kgf/cm^2
	I = 1.0 / 12;//���Ծ�cm^4
	A = 1 ;//��������cm^2
	L0 = 200 ;//��Ԫ0����200cm
	L1 = 200 * sqrt(2.0) ;//��Ԫ1����
	L2 = 200 * sqrt(2.0);//��Ԫ2����
	P = -10 ;//P=10 kgf
	//����洢�ռ�
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

	//���㵥��
	MatBeZero( ke0 , 6 , 6 ) ;//��������
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
	MatBeZero( ke1 , 6 , 6 ) ;//��������
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
	MatBeZero(ke2, 6, 6);//��������
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
	//1�ŵ�������ת��
	alpha = 5.0 / 4 * acos(-1.0) ;//����ת����=5/4*PI
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
	//2�ŵ�������ת��
	beta = 3.0 / 4 * acos(-1.0);//����ת����=3/4*PI
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
	//�鼯�ܸ�
	MatBeZero( K , 12 , 12 ) ;//�ܸ�������
	//0�ŵ����鼯
	//��Ԫ0����Ʊ
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
	//1�ŵ����鼯
	//��Ԫ1����Ʊ
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
	//2�ŵ����鼯
	//��Ԫ2����Ʊ
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

	//��ӡ�ܸ�
	cout<<"�ܸ�:"<<endl ;
	PrintMat( K , 12 , 12 ) ;
	//�غ�����
	MatBeZero( F , 12 , 1 ) ;
	F[ 4 ] = P ;
	//�߽���������
	//0�Žڵ�չ�(�Խ�Ԫ�ô�����)
	K[ 0 * 12 + 0 ] = 1e20 ;
	K[ 1 * 12 + 1 ] = 1e20 ;
	K[ 2 * 12 + 2 ] = 1e20 ;
	//2�Žڵ�չ�(�Խ�Ԫ�ô�����)
	K[ 6 * 12 + 6 ] = 1e20 ;
	K[ 7 * 12 + 7 ] = 1e20 ;
	K[ 8 * 12 + 8 ] = 1e20 ;
	//3�Žڵ�չ�(�Խ�Ԫ�ô�����)
	K[9 * 12 + 9] = 1e20;
	K[10 * 12 + 10] = 1e20;
	K[11 * 12 + 11] = 1e20;
	//���
	Gauss( K , V , F ,  12) ;
	//��ӡ�ڵ�λ������
	cout<<"�ڵ�λ������:"<<endl;
	PrintMat( V , 12 , 1 ) ;
	//���0�ŵ�Ԫ�˶���
	MatBeZero(F0, 6, 1);
	MatBeZero(V0, 6, 1);
	V0[0 * 1 + 0] = V[0 * 12 + 0]; V0[1 * 1 + 0] = V[1 * 1 + 0]; V0[2 * 1 + 0] = V[2 * 1 + 0];
	V0[3 * 1 + 0] = V[3 * 1 + 0]; V0[4 * 1 + 0] = V[4 * 1 + 0]; V0[5 * 1 + 0] = V[5 * 1 + 0];
	//F0=Ke0*V0
	MatMul(F0, ke0, V0, 6, 6, 1);
	cout << "0�ŵ�Ԫ�˶���:" << endl;
	PrintMat(F0, 6, 1);
	//���1�ŵ�Ԫ�˶���
	MatBeZero(F1, 6, 1);
	MatBeZero(V1, 6, 1);
	V1[0 * 1 + 0] = V[3* 12 + 0]; V1[1 * 1 + 0] = V[4 * 1 + 0]; V1[2 * 1 + 0] = V[5 * 1 + 0];
	V1[3 * 1 + 0] = V[6 * 1 + 0]; V1[4 * 1 + 0] = V[7 * 1 + 0]; V1[5 * 1 + 0] = V[8 * 1 + 0];
	//F1=Ke1*V1
	MatMul(F1, ke1, V1, 6, 6, 1);
	cout << "1�ŵ�Ԫ�˶���:" << endl;
	PrintMat(F1, 6, 1);
	//���2�ŵ�Ԫ�˶���
	MatBeZero(F2, 6, 1);
	MatBeZero(V2, 6, 1);
	V2[0 * 1 + 0] = V[3 * 12 + 0]; V2[1 * 1 + 0] = V[4 * 1 + 0]; V2[2 * 1 + 0] = V[5 * 1 + 0];
	V2[3 * 1 + 0] = V[9 * 1 + 0]; V2[4 * 1 + 0] = V[10* 1 + 0]; V2[5 * 1 + 0] = V[11 * 1 + 0];
	//F2=Ke2*V2
	MatMul(F2, ke2, V2, 6, 6, 1);
	cout << "2�ŵ�Ԫ�˶���:" << endl;
	PrintMat(F2, 6, 1);
	
	//�ͷ�ռ�ÿռ�
	delete []ke0;
	delete []ke1;
	delete []ke2;
	delete []K ;
	delete []F ;
	delete []V ;

	system("Pause");
	return 0;
}