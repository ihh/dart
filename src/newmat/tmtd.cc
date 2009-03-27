
//#define WANT_STREAM

#include "include.h"
#include "newmatap.h"

#include "tmt.h"

#ifdef use_namespace
using namespace NEWMAT;
#endif

ReturnMatrix Inverter(const CroutMatrix& X)
{
   Matrix Y = X.i();
   Y.Release();
   return Y.ForReturn();
}

void CircularShift(const Matrix& X1, int first, int last)
{
      Matrix X; UpperTriangularMatrix U1, U2;
      int n = X1.Ncols();

      // Try right circular shift of columns
      X = X1; QRZ(X, U1);
      RightCircularUpdateCholesky(U1, first, last);
      X = X1.Columns(1,first-1) | X1.Column(last)
         | X1.Columns(first,last-1) | X1.Columns(last+1,n);
      QRZ(X, U2);
      X = U1 - U2; Clean(X, 0.000000001); Print(X);

      // Try left circular shift of columns
      X = X1; QRZ(X, U1);
      LeftCircularUpdateCholesky(U1, first, last);
      X = X1.Columns(1,first-1) | X1.Columns(first+1,last)
         | X1.Column(first) | X1.Columns(last+1,n);
      QRZ(X, U2);
      X = U1 - U2; Clean(X, 0.000000001); Print(X);
}

class TestUpdateQRZ
{
   int m,n1,n2,n3;
   Matrix X1, X2, X3;
   MultWithCarry mwc;   // Uniform random number generator
public:
   void Reset();
   TestUpdateQRZ(int mx, int n1x, int n2x=0, int n3x=0)
      : m(mx), n1(n1x), n2(n2x), n3(n3x) { Reset(); }
   void DoTest();
   void ClearRow(int i)       { X1.Row(i) = 0.0; }
   void SetRow(int i, int j)  { X1.Row(i) = X1.Row(j); }
};



void trymatd()
{
   Tracer et("Thirteenth test of Matrix package");
   Tracer::PrintTrace();
   Matrix X(5,20);
   int i,j;
   for (j=1;j<=20;j++) X(1,j) = j+1;
   for (i=2;i<=5;i++) for (j=1;j<=20; j++) X(i,j) = (long)X(i-1,j) * j % 1001;
   SymmetricMatrix S; S << X * X.t();
   Matrix SM = X * X.t() - S;
   Print(SM);
   LowerTriangularMatrix L = Cholesky(S);
   Matrix Diff = L*L.t()-S; Clean(Diff, 0.000000001);
   Print(Diff);
   {
      Tracer et1("Stage 1");
      LowerTriangularMatrix L1(5);
      Matrix Xt = X.t(); Matrix Xt2 = Xt;
      QRZT(X,L1);
      Diff = L - L1; Clean(Diff,0.000000001); Print(Diff);
      UpperTriangularMatrix Ut(5);
      QRZ(Xt,Ut);
      Diff = L - Ut.t(); Clean(Diff,0.000000001); Print(Diff);
      Matrix Y(3,20);
      for (j=1;j<=20;j++) Y(1,j) = 22-j;
      for (i=2;i<=3;i++) for (j=1;j<=20; j++)
         Y(i,j) = (long)Y(i-1,j) * j % 101;
      Matrix Yt = Y.t(); Matrix M,Mt; Matrix Y2=Y;
      QRZT(X,Y,M); QRZ(Xt,Yt,Mt);
      Diff = Xt - X.t(); Clean(Diff,0.000000001); Print(Diff);
      Diff = Yt - Y.t(); Clean(Diff,0.000000001); Print(Diff);
      Diff = Mt - M.t(); Clean(Diff,0.000000001); Print(Diff);
      Diff = Y2 * Xt2 * S.i() - M * L.i();
      Clean(Diff,0.000000001); Print(Diff);
   }

   ColumnVector C1(5);
   {
      Tracer et1("Stage 2");
      X.ReSize(5,5);
      for (j=1;j<=5;j++) X(1,j) = j+1;
      for (i=2;i<=5;i++) for (j=1;j<=5; j++)
         X(i,j) = (long)X(i-1,j) * j % 1001;
      for (i=1;i<=5;i++) C1(i) = i*i;
      CroutMatrix A = X;
      ColumnVector C2 = A.i() * C1; C1 = X.i()  * C1;
      X = C1 - C2; Clean(X,0.000000001); Print(X);
   }

   {
      Tracer et1("Stage 3");
      X.ReSize(7,7);
      for (j=1;j<=7;j++) X(1,j) = j+1;
      for (i=2;i<=7;i++) for (j=1;j<=7; j++)
         X(i,j) = (long)X(i-1,j) * j % 1001;
      C1.ReSize(7);
      for (i=1;i<=7;i++) C1(i) = i*i;
      RowVector R1 = C1.t();
      Diff = R1 * X.i() - ( X.t().i() * R1.t() ).t(); Clean(Diff,0.000000001);
      Print(Diff);
   }

   {
      Tracer et1("Stage 4");
      X.ReSize(5,5);
      for (j=1;j<=5;j++) X(1,j) = j+1;
      for (i=2;i<=5;i++) for (j=1;j<=5; j++)
         X(i,j) = (long)X(i-1,j) * j % 1001;
      C1.ReSize(5);
      for (i=1;i<=5;i++) C1(i) = i*i;
      CroutMatrix A1 = X*X;
      ColumnVector C2 = A1.i() * C1; C1 = X.i()  * C1; C1 = X.i()  * C1;
      X = C1 - C2; Clean(X,0.000000001); Print(X);
   }


   {
      Tracer et1("Stage 5");
      int n = 40;
      SymmetricBandMatrix B(n,2); B = 0.0;
      for (i=1; i<=n; i++)
      {
         B(i,i) = 6;
         if (i<=n-1) B(i,i+1) = -4;
         if (i<=n-2) B(i,i+2) = 1;
      }
      B(1,1) = 5; B(n,n) = 5;
      SymmetricMatrix A = B;
      ColumnVector X(n);
      X(1) = 429;
      for (i=2;i<=n;i++) X(i) = (long)X(i-1) * 31 % 1001;
      X = X / 100000L;
      // the matrix B is rather ill-conditioned so the difficulty is getting
      // good agreement (we have chosen X very small) may not be surprising;
      // maximum element size in B.i() is around 1400
      ColumnVector Y1 = A.i() * X;
      LowerTriangularMatrix C1 = Cholesky(A);
      ColumnVector Y2 = C1.t().i() * (C1.i() * X) - Y1;
      Clean(Y2, 0.000000001); Print(Y2);
      UpperTriangularMatrix CU = C1.t().i();
      LowerTriangularMatrix CL = C1.i();
      Y2 = CU * (CL * X) - Y1;
      Clean(Y2, 0.000000001); Print(Y2);
      Y2 = B.i() * X - Y1; Clean(Y2, 0.000000001); Print(Y2);

      LowerBandMatrix C2 = Cholesky(B);
      Matrix M = C2 - C1; Clean(M, 0.000000001); Print(M);
      ColumnVector Y3 = C2.t().i() * (C2.i() * X) - Y1;
      Clean(Y3, 0.000000001); Print(Y3);
      CU = C1.t().i();
      CL = C1.i();
      Y3 = CU * (CL * X) - Y1;
      Clean(Y3, 0.000000001); Print(Y3);

      Y3 = B.i() * X - Y1; Clean(Y3, 0.000000001); Print(Y3);

      SymmetricMatrix AI = A.i();
      Y2 = AI*X - Y1; Clean(Y2, 0.000000001); Print(Y2);
      SymmetricMatrix BI = B.i();
      BandMatrix C = B; Matrix CI = C.i();
      M = A.i() - CI; Clean(M, 0.000000001); Print(M);
      M = B.i() - CI; Clean(M, 0.000000001); Print(M);
      M = AI-BI; Clean(M, 0.000000001); Print(M);
      M = AI-CI; Clean(M, 0.000000001); Print(M);

      M = A; AI << M; M = AI-A; Clean(M, 0.000000001); Print(M);
      C = B; BI << C; M = BI-B; Clean(M, 0.000000001); Print(M);
   }

   {
      Tracer et1("Stage 5");
      SymmetricMatrix A(4), B(4);
      A << 5
        << 1 << 4
        << 2 << 1 << 6
        << 1 << 0 << 1 << 7;
      B << 8
        << 1 << 5
        << 1 << 0 << 9
        << 2 << 1 << 0 << 6;
      LowerTriangularMatrix AB = Cholesky(A) * Cholesky(B);
      Matrix M = Cholesky(A) * B * Cholesky(A).t() - AB*AB.t();
      Clean(M, 0.000000001); Print(M);
      M = A * Cholesky(B); M = M * M.t() - A * B * A;
      Clean(M, 0.000000001); Print(M);
   }
   
   {
      Tracer et1("Stage 6");
      int N=49;
      int i;
      SymmetricBandMatrix S(N,1);
      Matrix B(N,N+1); B=0;
      for (i=1;i<=N;i++) { S(i,i)=1; B(i,i)=1; B(i,i+1)=-1; }
      for (i=1;i<N; i++) S(i,i+1)=-.5;
      DiagonalMatrix D(N+1); D = 1;
      B = B.t()*S.i()*B - (D-1.0/(N+1))*2.0;
      Clean(B, 0.000000001); Print(B);
   }
   {
      Tracer et1("Stage 7");
      // See if you can pass a CroutMatrix to a function
      Matrix A(4,4);
      A.Row(1) <<  3 <<  2 << -1 <<  4;
      A.Row(2) << -8 <<  7 <<  2 <<  0;
      A.Row(3) <<  2 << -2 <<  3 <<  1;
      A.Row(4) << -1 <<  5 <<  2 <<  2;
      CroutMatrix B = A;
      Matrix C = A * Inverter(B) - IdentityMatrix(4);
      Clean(C, 0.000000001); Print(C);
   }

   {
      Tracer et1("Stage 8");
      // Modification of Cholesky decomposition

      int i, j;

      // Build test matrix
      Matrix X(100, 10);
      MultWithCarry mwc;   // Uniform random number generator
      for (i = 1; i <= 100; ++i) for (j = 1; j <= 10; ++j)
         X(i, j) = 2.0 * (mwc.Next() - 0.5);
      Matrix X1 = X;     // save copy

      // Form sums of squares and products matrix and Cholesky decompose
      SymmetricMatrix A; A << X.t() * X;
      UpperTriangularMatrix U1 = Cholesky(A).t();

      // Do QR decomposition of X and check we get same triangular matrix
      UpperTriangularMatrix U2;
      QRZ(X, U2);
      Matrix Diff = U1 - U2; Clean(Diff, 0.000000001); Print(Diff);

      // Try adding new row to X and updating triangular matrix 
      RowVector NewRow(10);
      for (j = 1; j <= 10; ++j) NewRow(j) = 2.0 * (mwc.Next() - 0.5);
      UpdateCholesky(U2, NewRow);
      X = X1 & NewRow; QRZ(X, U1);
      Diff = U1 - U2; Clean(Diff, 0.000000001); Print(Diff);

      // Try removing two rows and updating triangular matrix
      DowndateCholesky(U2, X1.Row(20));
      DowndateCholesky(U2, X1.Row(35));
      X = X1.Rows(1,19) & X1.Rows(21,34) & X1.Rows(36,100) & NewRow; QRZ(X, U1);
      Diff = U1 - U2; Clean(Diff, 0.000000001); Print(Diff);

      // Circular shifts

      CircularShift(X, 3,6);
      CircularShift(X, 5,5);
      CircularShift(X, 4,5);
      CircularShift(X, 1,6);
      CircularShift(X, 6,10);
   }
   
   {
      Tracer et1("Stage 9");
      // Try updating QRZ, QRZT decomposition
      TestUpdateQRZ tuqrz1(10, 100, 50, 25); tuqrz1.DoTest();
      tuqrz1.Reset(); tuqrz1.ClearRow(1); tuqrz1.DoTest();
      tuqrz1.Reset(); tuqrz1.ClearRow(1); tuqrz1.ClearRow(2); tuqrz1.DoTest();
      tuqrz1.Reset(); tuqrz1.ClearRow(5); tuqrz1.ClearRow(6); tuqrz1.DoTest();
      tuqrz1.Reset(); tuqrz1.ClearRow(10); tuqrz1.DoTest();
      TestUpdateQRZ tuqrz2(15, 100, 0, 0); tuqrz2.DoTest();
      tuqrz2.Reset(); tuqrz2.ClearRow(1); tuqrz2.DoTest();
      tuqrz2.Reset(); tuqrz2.ClearRow(1); tuqrz2.ClearRow(2); tuqrz2.DoTest();
      tuqrz2.Reset(); tuqrz2.ClearRow(5); tuqrz2.ClearRow(6); tuqrz2.DoTest();
      tuqrz2.Reset(); tuqrz2.ClearRow(15); tuqrz2.DoTest();
      TestUpdateQRZ tuqrz3(5, 0, 10, 0); tuqrz3.DoTest();
      
   }
   
//   cout << "\nEnd of Thirteenth test\n";
}

void TestUpdateQRZ::Reset()
{
   int i, j;
   // set up Matrices
   X1.ReSize(m, n1); X2.ReSize(m, n2); X3.ReSize(m, n3);
   for (i = 1; i <= m; ++i) for (j = 1; j <= n1; ++j)
      X1(i, j) = 2.0 * (mwc.Next() - 0.5);
   for (i = 1; i <= m; ++i) for (j = 1; j <= n2; ++j)
      X2(i, j) = 2.0 * (mwc.Next() - 0.5);
   for (i = 1; i <= m; ++i) for (j = 1; j <= n3; ++j)
      X3(i, j) = 2.0 * (mwc.Next() - 0.5);
}

void TestUpdateQRZ::DoTest()
{
   Matrix XT1 = X1.t(), XT2 = X2.t(), XT3 = X3.t();
   Matrix X = X1 | X2 | X3;
   SymmetricMatrix SM; SM << X * X.t();
   LowerTriangularMatrix L, LX, LY, L0;
   QRZT(X, L);
   X = X1 | X2 | X3;
   LY.ReSize(m); LY = 0.0;
   UpdateQRZT(X, LY);
   QRZT(X1, LX); UpdateQRZT(X2, LX); UpdateQRZT(X3, LX);
   Matrix Diff = L * L.t() - SM; Clean(Diff, 0.000000001); Print(Diff);
   Diff = SM - LX * LX.t(); Clean(Diff, 0.000000001); Print(Diff);
   Diff = SM - LY * LY.t(); Clean(Diff, 0.000000001); Print(Diff);
   UpperTriangularMatrix U, UX, UY;
   X = XT1 & XT2 & XT3;
   QRZ(X, U);
   Diff = U.t() * U - SM; Clean(Diff, 0.000000001); Print(Diff);
   UY.ReSize(m); UY = 0.0;
   X = XT1 & XT2 & XT3;
   UpdateQRZ(X, UY);
   Diff = SM - UY.t() * UY; Clean(Diff, 0.000000001); Print(Diff);
   QRZ(XT1, UX); UpdateQRZ(XT2, UX); UpdateQRZ(XT3, UX);
   Diff = SM - UX.t() * UX; Clean(Diff, 0.000000001); Print(Diff);
}   
