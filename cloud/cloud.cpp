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
            p[i*m+i] = 1.0;
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

void shape(matrix A)
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
    //cout<<"old_key shape: ("<<old_key.m<<" "<<old_key.n<<")"<<endl;
    //cout<<"T "; shape(T);
    //cout<<"max of new S: ";
    //cout<<max(T)<<"*****************"<<endl;
    return keySwitchMatrix(old_key,T);
}

/*return c:
This is for sever side to compute the inner product on ciphers c1 and c2, M is the
transform matrix to calculate the inner product.

@M: Transform matrix generated by client to
*/
const matrix innerProduct(const matrix& M, const matrix& c1, const matrix& c2,double w)
{
//    cout<<"M shape: ("<<M.m<<" "<<M.n<<")"<<endl;
    matrix tmp = (vectorize(c1*transpose(c2))/w);
//    cout<<"c1*c2/w shape: ("<<tmp.m<<" "<<tmp.n<<")"<<endl;
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
			out<<","<<a.p[i*a.n+j];
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
		a.p[i*n+j] = atof(row[j].c_str());
	}
	return a;
}
int main(int argc, char * argv[]){

    if(argc<1){cout<<"args < 2 ..."<<endl; return 0;}
    //double w = 1<<16;
    double w = (1<<30)*(1.0)*(1<<2);

    //int n = atoi(argv[1]);
	
	/*
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
	matrix c2 = encrypt(T2,b,w);
	int t = clock();
	cout<<"Encryption time: "<<double(t-s)/CLOCKS_PER_SEC/2<<"s\n";
	matrix S1= getSecretKey(T1), S2=getSecretKey(T2);
	*/
	if(!strcmp(argv[1],"T_alice")){
		matrix T(1,extension);
		//matrix S = getSecretKey(T);
		save(T,"T_alice.csv");
	}
	/*
	s = clock();
	matrix M = innerProductClient(S1,S2,matrix(n,'i'),T);
	t = clock();
	cout<<"M generation time: "<<double(t-s)/CLOCKS_PER_SEC<<endl;
    matrix H(n,'i');
	*/
	if(!strcmp(argv[1],"Product"))
	{
		matrix c1 = load("../c1.csv");
		cout<<"c1: "<<endl;
		c1.display();
		matrix c2 = load("../c2.csv");
		cout<<"c2: "<<endl;
		c2.display();
		matrix M = load("../M.csv");
		
		int s = clock();
    matrix innerProd = innerProduct(M,c1,c2,w);	
	cout<<"InnerProd cipher: "<<endl;
	innerProd.display();
    int t = clock();
    cout<<"Inner Product Calculation time: "<<double(t-s)/CLOCKS_PER_SEC<<endl;	
	matrix S = getSecretKey(load("T_alice.csv"));
    matrix dec = decrypt(S,innerProd,w);
	cout<<"Inner Product:\n";
	dec.display();	
	}

	if(!strcmp(argv[1],"3"))
	{
		/*length of plaintext*/
		int n = atoi(argv[2]);
		extension = 10;
		matrix T_cloud = load("T_alice.csv");
		matrix Snew = getSecretKey(T_cloud);
		matrix M1 = load("../alice/M1.csv");
		//cout<<"M1 shape:";
		//shape(M1);
		matrix M2 = load("../bob/M2.csv");
		//cout<<"M2 shape:";
		//shape(M2);
		/*Construct M3*/
		matrix I(n,'i');
		matrix up_right(n,extension);
		zero(up_right);
		matrix up = A_concatnate_B(I,up_right);

		matrix down_left(extension,n);
		zero(down_left);
		
		matrix T1 = load("../alice/T1.csv");
		matrix T2 = load("../bob/T2.csv");
		matrix tt = transpose(T1)*T2; 
		//cout<<"T1*T2 shape: "; shape(tt);
		//cout<<"down_left shape: "; shape(down_left);
		matrix S3 = A_on_B(up,A_concatnate_B(down_left,transpose(T1)*T2));
		//cout<<"S3 is: "<<endl;
		//S3.display();
		int s = clock();
		matrix M3 = keySwitchMatrix(transpose(vectorize(S3)),T_cloud);
		/*
		cout<<"M3 shape:";
		shape(M3);
		cout<<"S*M =: "<<endl;
		cout<<"Snew shape:";
		shape(Snew);
		cout<<"S*M shape:";
		*/
		matrix M = M1+M2+M3;
		int t = clock();
		cout<<setprecision(2)<<double(t-s)/CLOCKS_PER_SEC<<endl;
		/*
		matrix SM = Snew*M;
		shape(SM);
		cout<<"S1*S2:";
		
		matrix S1 = getSecretKey(T1);
		matrix S2 = getSecretKey(T2);
		shape(transpose(S1)*S2);
		matrix S = transpose(vectorize(transpose(S1)*S2));
		matrix Sb = getBitMatrix(S);
		cout<<"Sb shape: ";
		shape(Sb);
		/*
		if(SM == Sb)
			cout<<"hahahahaha"<<endl;
		else
			cout<<"wocao"<<endl;
		*/
		matrix c1 = load("../alice/c1.csv");
		matrix* cipher_matrix1 = new matrix[c1.n];
		for(int j=0;j<c1.n;j++)
			cipher_matrix1[j] = matrix(c1.m,1);
		for(int i=0;i<c1.m;i++)
			for(int j=0;j<c1.n;j++)
				{
					cipher_matrix1[j] = matrix(c1.m,1);
					cipher_matrix1[j].p[i] = c1.p[i*c1.n+j];
				}
		int cnt_a = c1.n;

		matrix c2 = load("../bob/c2.csv");
		matrix* cipher_matrix2 = new matrix[c2.n];

		for(int j=0;j<c2.n;j++)
			cipher_matrix2[j] = matrix(c2.m,1);
		for(int i=0;i<c2.m;i++)
			for(int j=0;j<c2.n;j++)
			{
				cipher_matrix2[j].p[i] = c2.p[i*c2.n+j];
			}
		int cnt_b = c2.n;
		
		//int cnt =0;
		int s1;
		double cnt1=0,cnt2=0,cnt3=0;
		// cal interweaved inner product
		//s = clock();
		//double kernel_inter[cnt_a][cnt_b];
		for(int i=0;i<cnt_a;i++)
			for(int j=0;j<cnt_b;j++)
			{
			//	cnt++;
			s1 = clock();
				matrix cp = vectorize(cipher_matrix1[i]*transpose(cipher_matrix2[j]))/w;
			cnt1 += double(clock()-s1)/CLOCKS_PER_SEC;
			s1 = clock();
				matrix newc = keySwitch(M,cp);
			cnt2 += double(clock()-s1)/CLOCKS_PER_SEC;
			s1 = clock();
				matrix product = decrypt(Snew,newc,w);
			cnt3 += double(clock()-s1)/CLOCKS_PER_SEC;
		/*
		cout<<"product is:";
		product.display();
		*/
		//		kernel_inter[i][j] = product.p[0]; 
			}
		//t = clock();
		//cout<<setprecision(2)<<double(t-s)/CLOCKS_PER_SEC<<endl;
		//cout<<"M: "<<M.m<<" "<<M.n<<endl;
		/**
		 * cal inner products of alice's data
		 */
		matrix Ma = load("../alice/Ma.csv");
		
		//int cnt_aa=0;
		//s = clock();
		//cout<<"Ma: "<<Ma.m<<" "<<Ma.n<<endl;
		double cnta1=0,cnta2=0,cnta3=0;
		//double kernel_alice[cnt_a][cnt_a];
		for(int i=0;i<cnt_a;i++)
		{
			for(int j=i;j<cnt_a;j++)
			{
		//	cnt_aa++;
			s1 = clock();
			matrix cp = vectorize(cipher_matrix1[i]*transpose(cipher_matrix1[j]))/w;
			cnta1 += double(clock()-s1)/CLOCKS_PER_SEC;
			s1 = clock();
			matrix newc = keySwitch(Ma,cp);
			cnta2 += double(clock()-s1)/CLOCKS_PER_SEC;
			s1 = clock();
			matrix product = decrypt(Snew,newc,w);
		//	kernel_alice[i][j] = product.p[0];
			cnta3 += double(clock()-s1)/CLOCKS_PER_SEC;
			}
		}
		//t = clock();
		//cout<<setprecision(2)<<double(t-s)/CLOCKS_PER_SEC<<endl;
		
		int cnt_bb=0;
		//call inner products of bob's data
		matrix Mb = load("../bob/Mb.csv");
		//s = clock();
		double cntb1=0,cntb2=0,cntb3=0;
		//double kernel_bob[cnt_b][cnt_b];
		for(int i=0;i<cnt_b;i++)
		{
			for(int j=i;j<cnt_b;j++)
			{
			//cnt_bb++;
			s1 = clock();
			matrix cp = vectorize(cipher_matrix2[i]*transpose(cipher_matrix2[j]))/w;
			cntb1 += double(clock()-s1)/CLOCKS_PER_SEC;
			s1 = clock();
			matrix newc = keySwitch(Mb,cp);
			cntb2 += double(clock()-s1)/CLOCKS_PER_SEC;
			s1 = clock();
			matrix product = decrypt(Snew,newc,w);
		//	kernel_bob[i][j] = product.p[0];
			cntb3 += double(clock()-s1)/CLOCKS_PER_SEC;
			}
		}
		cout<<setprecision(2)<<cnt1+cnt2<<endl;
		cout<<setprecision(2)<<cnt3<<endl;
		cout<<setprecision(2)<<cnta1+cnta2<<endl;
		cout<<setprecision(2)<<cnta3<<endl;
		cout<<setprecision(2)<<cntb1+cntb2<<endl;
		cout<<setprecision(2)<<cntb3<<endl;
		//t = clock();
		//cout<<setprecision(2)<<double(t-s)/CLOCKS_PER_SEC<<endl;
		//cout<<cnt<<endl;
		//cout<<cnt_aa<<endl;
		//cout<<cnt_bb<<endl;
	}
	return 0;

    /*Weighted Inner Products*/
    /* Setting w >> |S|,where S is any secret_key (after encryption or keySwitch)
        to make sure decryption can be successful.
    */
    /*
    int s = clock();
    matrix c = encrypt(T,a,w);
    int t = clock();
    cout<<double(t-s)/CLOCKS_PER_SEC<<endl;
    */
    
//  for inner product , row number of T should equals to S*, which is always 1

    //matrix old_key = transpose(vectorize(transpose(S1)*H*S2));
    //matrix innerProduct= decrypt(S,innerProd,w);
    //innerProduct.display();
    /*
    matrix cipher = encrypt(T,plain,w);
    cout<<"enced: "<<endl;
    cipher.display();

    matrix dec = decrypt(getSecretKey(T),cipher,w);
    cout<<"deced: "<<endl;
    dec.display();
    */
    return 0;
}
