//
// Copyright (c) 2015 by Xin Wang.
// All Rights Reserved.
//


#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h> /*srand, rand*/
#include <time.h>   /*time*/
#include <assert.h>
#include <string.h>
#include <vector>
#include <iomanip>
#include "mkl.h"

using namespace std;

/**
    matrix class

    This is a matrix class, mainly supporting fast matrix multiplication by call the
    cblas_dgemm function in mkl.
 */

class matrix{
public:
    /**
        member variables
        @m: number of row in the matrix
        @n: number of column in the matrix
        @p: A pointer to a 1-D double array to represent 2-D matrix, with lenth m*n*sizeof(double)
     */
    int m,n;
    double *p;
    //matrix(){}
    /**
        Constructor with default args.
        init a matrix with shape(row,col), elements of it are randomized.
        @row: row number of matrix
        @col: col number of matrix
     */
    matrix(int row = 2, int col = 2)
    {
        m= row;
        n = col;
        p = (double*)mkl_malloc(m*n*sizeof(double),64);
        srand(time(NULL));
        for(int i=0;i<m;i++)
        {
            for(int j=0;j<n;j++){
                /* An integer value between 0 and RAND_MAX.
                  RAND_MAX is a constant defined in <cstdlib>*/
                p[i*n+j] = rand()%(1<<11);
            }
        }
    }
    /**
        Copy Constructor
     */
    matrix(const matrix& a) //reference must
    {
        //matrix(a.m,a.n);
        /*First initialize "this" object(memory allocation),
        or it will lead to undefined behavior(bus error:10)*/
        m= a.m;
        n = a.n;
        p = (double*)mkl_malloc(m*n*sizeof(double),64);

        for(int i=0;i<m;i++)
            for(int j=0;j<n;j++)
                p[i*n+j] = a.p[i*n+j];
    }
    /**
        Identity matrix Constructor
        @d: dimension of the Identity
        @c: char should be set to 'i', tell the Constructor to generate identity matrix
     */
    matrix(int d, char c)
    {
        assert(c == 'i');
        m = n = d;
        p = (double*)mkl_malloc(m*n*sizeof(double),64);
        memset(p,0,m*n*sizeof(double));
        for(int i=0;i<m;i++)
		{
				p[i*m+i] = 1.0;
		}
    }
    /**
        Destructor
        Call mkl_free to free array pointer p.
     */
    ~matrix()
    {
        mkl_free(p);
    }
    /**
        Transpose a matrix
     */
    friend matrix transpose (const matrix&);
    /**
        Overloading operator * between matries
        @args: two matries as args, number of cols of the first matrix equals to number
        rows of the second matrix.
     */
    friend matrix operator * (const matrix&, const matrix&);
    /*Overloading operator * between matrix and double
     */
    friend matrix operator * (const matrix&, const double);
    /*Overloading operator * between double and matrix
     */
    friend matrix operator * (const double , const matrix&);
    /*Overloading operator /
     */
    const matrix& operator / (const double );
    friend matrix operator + (const matrix&, const matrix&);
    friend matrix operator - (const matrix&, const matrix&);

    friend matrix A_on_B (const matrix&, const matrix&);
    friend matrix A_concatnate_B (const matrix&, const matrix&);
    friend bool operator == (const matrix& A, const matrix& B);
    /*display
      Display the matrix
    */
    void display()
    {
//        cout<<"The matrix is";
        for(int i=0;i<m;i++){
            //cout<<endl;
            for(int j=0;j<n;j++){
                cout<<long (p[i*n+j])<<" ";
            }
			cout<<endl;
        }
        //cout<<endl;
    }

};

void shape(matrix & A)
{
    cout<<"shape is: ("<<A.m<<" "<<A.n<<")"<<endl;
}

void zero(matrix & A)
{
	memset(A.p,0,A.m*A.n*sizeof(double));
}
matrix transpose(const matrix& A)
{
    matrix res(A.n,A.m);
    for(int i=0;i<A.m;i++)
        for(int j=0;j<A.n;j++)
            res.p[j*A.m+i] = A.p[i*A.n+j];
    return res;
}

matrix A_on_B (const matrix& A, const matrix& B)
{
    assert(A.n == B.n);
    int m = A.m+B.m, n = A.n;
    matrix res(m,n);
    for(int i=0;i<A.m;i++)
        for(int j=0;j<n;j++)
            res.p[i*n+j] = A.p[i*n+j];
    for(int i=0;i<B.m;i++)
        for(int j=0;j<n;j++)
            res.p[(i+A.m)*n+j] = B.p[i*n+j];
    return res;
}

matrix A_concatnate_B (const matrix& A, const matrix& B)
{
    assert(A.m == B.m);
    int m = A.m,n=A.n+B.n;
    matrix *res = new matrix(m,n);
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++){
            if(j<A.n)  res->p[i*n+j] = A.p[i*A.n+j];
            else  res->p[i*n+j] = B.p[i*B.n+j-A.n];
        }
    return *res;
}

bool operator == (const matrix& A, const matrix& B)
{
    assert(A.m==B.m && A.n == B.n);
    for(int i=0;i<A.m;i++)
        for(int j=0;j<A.n;j++)
            if(A.p[i*A.n+j] != B.p[i*A.n+j]) return false;
    return true;
}

matrix operator + (const matrix& A , const matrix& B)
{
    assert(A.m==B.m && A.n==B.n);
    int m=A.m, n=A.n;
    matrix res(m,n);
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            res.p[i*n+j] = A.p[i*n+j] + B.p[i*n+j];
    return res;
}

matrix operator - (const matrix& A , const matrix& B)
{
    int m=A.m,n=A.n;
    matrix res = matrix(m,n);
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            res.p[i*n+j] = A.p[i*n+j] - B.p[i*n+j];
    return res;
}

matrix operator * (const matrix & A , const matrix & B)
{
    assert(A.n == B.m);
    matrix * res= new matrix (A.m,B.n);
    int m,n,p;
    m = A.m;
    p = A.n;
    n = B.n;
    double alpha = 1.0, beta = 0.0;
    //int a = clock();
    //double s_initial,s_elapsed;
    //s_initial = dsecnd();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,m,n,p,alpha,A.p,p,B.p,n,beta,(*res).p,n);
    //int b = clock();
    //cout<<double(b-a)/CLOCKS_PER_SEC<<endl;
    return *res;
}


matrix operator *(const matrix& A,const double a)
{
    matrix res(A);
    for(int i=0;i<A.m;i++)
        for(int j=0;j<A.n;j++)
            res.p[i*A.n+j] =A.p[i*A.n+j]*a;
    return res;
}

matrix operator *(const double a,const matrix& A)
{
    return operator*(A,a);
}

const matrix& matrix::operator / (double div)
{
    int m = this->m, n = this->n;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            this->p[i*n+j]/=div;
    return *this;
}

int l = 60;
int extension = 10;

double max(const matrix& A)
{
    double tmp = -(1<<20);
    for(int i=0;i<A.m*A.n;i++)
        if(A.p[i]>tmp) tmp=A.p[i];
    return tmp;
}

/**
    get c*

    This method gets the binary represented matrix of c by representing every element
    in c with l bit( 2^l>|c|,where |c| is the max abs value of c, so that no loss
    *happens ). For every element in c*, possible values are 1, 0, -1.

    Examples:

    Represent 7 and -7 with l=4 bits:
    7  --> [1,1,1,0]
    -7 --> [-1,-1,-1,0]

    @c: ciphertext
 */
matrix getBitVector(matrix c)
{
    int n = c.m;
    assert(c.n==1);
    matrix *res = new matrix(n*l,1);
    memset(res->p,0,n*l*sizeof(double));
    for(int i=0;i<n;++i){  //1 element to l-bit
        int sign = 1;
        if(c.p[i]<0) sign = -1;
        for(int j=0;j<l;++j){
            res->p[i*l+j] = sign * (int(c.p[i])&1);
            if(c.p[i])
            c.p[i]/=2;
        }
    }
    return *res;
}

/**
    get S*

    This method gets the binary represented matrix of secret_key by representing every
    elementin secret_key with l bit( 2^l>|secret_key|,where |secret_key| is the max
    abs value of secret_key, so that no loss happens. ). For every element in c*,
    possible values are 1, 0, -1.

    Examples:

    Represent 7 and -7 with l=4 bits:
    7  --> [1,1,1,0]
    -7 --> [-1,-1,-1,0]

    @secret_key: of its literal meaning
 */
matrix getBitMatrix(const matrix& secret_key)
{
    int m = secret_key.m , n = secret_key.n;
    matrix * res = new matrix(m,n*l);
    memset(res->p,0,m*n*l*sizeof(double));             //memset to 0
    for(int i=0;i<m;++i)
        for(int j=0;j<n;++j)
        {
            long tmp = secret_key.p[i*n+j];         //tmp should be long
            for(int k=0;k<l;++k)
                {
                    if(tmp==0) break;
                    res->p[i*n*l+j*l+k] = tmp<<k;
                }
        }
    return *res;
}

/**
    get M
    S: original secret_key
    T: T in new secret_key = [I,T]
 */
const matrix keySwitchMatrix(const matrix& S,const matrix& T)
{
    // bit string of secret_key
    matrix bitwise_secret_key = getBitMatrix(S);
    //cout<<"S*"; shape(bitwise_secret_key);
    // get random matrix A
    matrix A(T.n,bitwise_secret_key.n);
    return A_on_B(bitwise_secret_key-(T*A),A);
}

/**
    return c = M * c*;
 */
const matrix keySwitch(const matrix & M, const matrix & c)
{
		return M * getBitVector(c);
}

/**
    S = [I,T]
 */
const matrix encrypt(const matrix& T,const matrix& x,const double w)
{
    matrix c = keySwitch(keySwitchMatrix(matrix(T.m,'i'),T),w * x);
    return c;
}

/**
    return the nearest integer in the form of double.
 */
double nearestInt(double a)
{
    if((a-int(a))>=0.5000) return 1.0+int(a);
    else return double(int(a));
}
/**/
const matrix decrypt(const matrix & S, const matrix & c, const double w)
{
    matrix Sc = S*c;
    for(int i=0;i<Sc.m;i++)
        for(int j=0;j<Sc.n;j++)
            Sc.p[i*Sc.n+j] = nearestInt(Sc.p[i*Sc.n+j]/w);
    return Sc;
}

const matrix getSecretKey(const matrix& T)
{
    return A_concatnate_B(matrix(T.m,'i'),T);
}

/*  return M corresponding to Transform G
    send M to server to transform
    S: current secret_key
    G: Transform matrix
    T: T in new secret_key(generated by client)
*/
const matrix linearTransformClient(const matrix &T, const matrix& G,const matrix S)
{
    return keySwitchMatrix(G*S,T);
}

/*
    server-side transform
*/
const matrix linearTransform(const matrix& M, const matrix& c)
{
    return M*getBitVector(c);
}
/*
*/
matrix vectorize(const matrix& H)
{
    int m=H.m,n=H.n;
    matrix res(m*n,1);
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            res.p[j*m+i]= H.p[i*n+j];
    return res;
}
/*  return M
    S1: secret_key1
    S2: secret_key2
    H: weighted matrix
    T: T in new secret_key
*/
const matrix innerProductClient(const matrix& S1, const matrix& S2,const matrix& H, matrix& T)
{
    matrix old_key = transpose(vectorize(transpose(S1)*H*S2));
    return keySwitchMatrix(old_key,T);
}

/*return c:
This is for sever side to compute the inner product on ciphers c1 and c2, M is the
transform matrix to calculate the inner product.

@M: Transform matrix generated by client to
*/
const matrix innerProduct(const matrix& M, const matrix& c1, const matrix& c2,double w)
{
    matrix tmp = (vectorize(c1*transpose(c2))/w);
    return M*getBitVector(tmp);
}

void usage()
{
	cout<<"Usage: ./project [dimension] [opration]"<<endl;
}
/**
	save a matrix in a file in csv format, with delimiter "," 
 */
void save(matrix a,char * filename)
{
	ofstream out;
	out.open(filename);
	out<<fixed<<setprecision(0);
	out<<a.m<<","<<a.n<<endl;
	for(int i=0;i<a.m;i++)
	{
		out<<a.p[i*a.n+0];
		for(int j=1;j<a.n;j++)
		{
			out<<","<<a.p[i*a.n+j];
	//		cout<<","<<a.p[i*a.n+j];
		}
		out<<endl;
	}
	out.close();
}


class CSVRow
{
    public:
        string const& operator[](size_t index) const
        {
            return m_data[index];
        }
        size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(istream& str)
        {
            string         line;
            getline(str,line);

            stringstream   lineStream(line);
            string         cell;

            m_data.clear();
            while(getline(lineStream,cell,','))
            {
                m_data.push_back(cell);
            }
        }
    private:
        vector<string>    m_data;
};

istream& operator>> (istream& str, CSVRow& data)
{
	data.readNextRow(str);
	return str;
}
/**
	load the saved matrix in a file
 */
matrix load(char * filename)
{
	ifstream in(filename);
	CSVRow row;
	in>>row;
	int m=atoi(row[0].c_str()), n = atoi(row[1].c_str());
	matrix a(m,n);
	for(int i=0;i<m;i++)
	{
		in>>row;
		for(int j=0;j<n;j++)
		a.p[i*n+j] = atoi(row[j].c_str());
	}
	return a;
}
int main(int argc, char * argv[]){
    if(argc<1){cout<<"args < 2 ..."<<endl; return 0;}
	
    double w = (1<<30)*(1.0)*(1<<2);
	
    //double w = 1<<16;
    int n = atoi(argv[1]);
    int cnt = atoi(argv[2]);
	/*	
	matrix a(n,1);
	cout<<"a is:"<<endl;
	a.display();
	matrix c = getBitVector(a);
	cout<<"bits is:"<<endl;
	c.display();
	*/
/*
	extension=2;
		
	matrix a(n,1);
	matrix b(n,1);
	cout<<"vector a:\n";
	transpose(a).display();
	cout<<"vector b:\n";
	transpose(b).display();
	cout<<"Plain Dot Product is: \n";
	(transpose(a)*b).display();
	matrix T1(n,extension),T2(n,extension);
	int s = clock();

	matrix c1 = encrypt(T1,a,w);
	cout<<"c1: "<<endl;
	c1.display();

	save(c1,"c1.csv");
	matrix c2 = encrypt(T2,b,w);
	save(c2,"c2.csv");

	int t = clock();
	cout<<"Encryption time: "<<double(t-s)/CLOCKS_PER_SEC/2<<"s\n";
	matrix S1= getSecretKey(T1), S2=getSecretKey(T2);
	//matrix T(1,extension);
	matrix T = load("../cloud/T_alice.csv");
	matrix S = getSecretKey(T);
	s = clock();
	matrix M = innerProductClient(S1,S2,matrix(n,'i'),T);
	save(M,"M.csv");
	
	t = clock();
	cout<<"M generation time: "<<double(t-s)/CLOCKS_PER_SEC<<endl;
    matrix H(n,'i');
    s = clock();
    matrix innerProd = innerProduct(M,c1,c2,w);	
    t = clock();
	cout<<"innerProd cipher:"<<endl;
	innerProd.display();
    cout<<"Inner Product Calculation time: "<<double(t-s)/CLOCKS_PER_SEC<<endl;	
    matrix dec = decrypt(getSecretKey(T),innerProd,w);
	cout<<"Inner Product:\n";
	dec.display();
	//*/

	matrix* data = new matrix[cnt];	
	for(int i=0;i<cnt;i++)
		data[i] = matrix(n,1);
	//cout<<"a is:"<<endl;
	//a.display();
	//save(a,"p1.csv");
	matrix T(n,extension);
	
	matrix *cipher = new matrix[cnt];
	int s = clock();
	for(int i=0;i<cnt;i++)
		cipher[i] = encrypt(T,data[i],w);
	int t = clock();
    cout<<setprecision(2)<<double(t-s)/CLOCKS_PER_SEC<<endl;	
	//cout<<"c1"<<endl;
	//c1.display();
	// save cipher matrix
	
	matrix cipher_matrix(n+extension,cnt);
	for(int i=0;i<n+extension;i++)
		for(int j=0;j<cnt;j++)
			cipher_matrix.p[i*cnt+j] = cipher[j].p[i];

	save(cipher_matrix,"c1.csv");

	matrix S = getSecretKey(T);
	matrix M1_tmp = transpose(S);	
	zero(M1_tmp);
	matrix r(n+extension,extension);
	zero(r);
	
	matrix S1 = A_concatnate_B(M1_tmp,r);
	matrix T_cloud = load("../cloud/T_alice.csv");
	// cal key switch matrix M1
	s = clock();
	matrix M1 = keySwitchMatrix(transpose(vectorize(S1)),T_cloud);
	t = clock();
	cout<<setprecision(2)<<double(t-s)/CLOCKS_PER_SEC<<endl;
	// save M1
	save(M1,"M1.csv");

	// cal key switch matrix M for alice's data
	s = clock();
	matrix M = innerProductClient(S,S,matrix(n,'i'),T_cloud);
	t = clock();
	cout<<setprecision(2)<<double(t-s)/CLOCKS_PER_SEC<<endl;
	// save M
	save(M,"Ma.csv");
	//random matrix R
	matrix R(T.m,T.n); 
	save(T+R,"T1.csv");
	save(R,"R1.csv");
	
	return 0;

    /*Weighted Inner Products*/
    /* Setting w >> |S|,where S is any secret_key (after encryption or keySwitch)
        to make sure decryption can be successful.
    */
    
//  for inner product , row number of T should equals to S*, which is always 1

    return 0;
}
