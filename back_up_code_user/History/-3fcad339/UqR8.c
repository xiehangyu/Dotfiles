    #include "../include/MatrixTool3.h"


void Debug_ShowMatrix(int n, int m, double **M, int exit_)
{
    // for debug
    for (int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            printf("%f   ", M[i][j]);
        }
        printf("\n");
    }
    if (1 == exit_){
        exit(EXIT_SUCCESS);
    }
}
    

void Debug_ShowVector(int n, double *Vect, int exit_)
{
    // for debug
    for (int i = 0; i < n; i++)
    {
        printf("%f   ", Vect[i]);
    }
    printf("\n");
    if (1 == exit_){
        exit(EXIT_SUCCESS);
    }
}


double** MatrixGenerator(int n, int m, void(*Init)(int, int, double**))
{
    // detect input error
    if (n<=0 || m<=0){
        printf("Matrix parameters input error!!!");
        exit(EXIT_FAILURE);
    }

    // generate Matrix with malloc
    double **Matrix = (double **)malloc(sizeof(double *)*n);
    for (int i=0;i<n;i++){
        Matrix[i] = (double *)malloc(sizeof(double)*m);
    }

    Init(n, m, Matrix);

    return Matrix;
}


void FreeMatrixMemory(double **Matrix, int n)
{
    for (int i=0;i<n;i++){
        free(Matrix[i]);
    }
    free(Matrix);
}


void MatrixInit_Rand_0_1(int n, int m, double **M)
{
    // set random number seed
    srand((unsigned) time(NULL));

    for (int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            M[i][j] = (double)rand() / (double)RAND_MAX;
        }
    }
}


void MatrixInit_0(int n, int m, double **M)
{
    for (int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            M[i][j] = 0.;
        }
    }
}


void VectInit_Same(int n, double *Vector, double value)
{
    for (int i=0;i<n;i++){
        Vector[i] = value;
    }
}


double Distance(int n, double *X1, double *X0)
{
    double distance = 0.;

    for (int i=0;i<n;i++){
        distance = distance + (X1[i] - X0[i])*(X1[i] - X0[i]);
    }
    distance = sqrt(distance);

    return distance;
}


void VectNorm(int n, double *Vector)
{
    double VectorLenth = 0.; // |Vector|

    for (int i=0;i<n;i++){
        VectorLenth = VectorLenth + Vector[i]*Vector[i];
    }
    VectorLenth = sqrt(VectorLenth);

    for (int i=0;i<n;i++){
        Vector[i] = Vector[i] / VectorLenth;
    }  
}


void MatrixMutiply(int n, int m, double **Matrix1, double **Matirx2, double **Matrix3)
{
    // buffer
    double **Matrix12 = MatrixGenerator(n, n, MatrixInit_0);

    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            for (int k=0;k<m;k++){
                Matrix12[i][j] = Matrix12[i][j] + Matrix1[i][k]*Matirx2[k][j];
            }
        }
    }

    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            Matrix3[i][j] = Matrix12[i][j];
        }
    }    

    FreeMatrixMemory(Matrix12, n);
}


void MatrixTranspose(double n, double m, double **Matrix, double **MatrixT)
{
    for (int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            MatrixT[j][i] = Matrix[i][j];
        }
    }
}


void GaussElim_PartialPivoting(int n, double **A, double *X, double *b, int DetIs0)
{
    /*suppose A as nxn here*/

    // 上三角化
    for (int i=0; i<n; i++)
    {
        // Select column pivot
        int iMax = i;
        for (int k=i+1;k<n;k++){
            if (fabs(A[iMax][i]) < fabs(A[k][i])){
                iMax = k;
            }
        }
        for (int j=i; j<n; j++)
        {
            double term;
            // Permuting(A[i][j], A[iMax][j]);
            term = A[i][j];
            A[i][j] = A[iMax][j];
            A[iMax][j] = term;            
        }
        // Permuting(b[i], b[iMax]);
        double term = b[i];
        b[i] = b[iMax];
        b[iMax] = term;
        
        // Elimination.
        for (int k=i+1; k<n; k++){
            A[k][i] = A[k][i]/A[i][i];
            for (int j=i+1; j<n; j++){
                A[k][j] = A[k][j] - A[k][i]*A[i][j];
            }
            b[k] = b[k] - A[k][i]*b[i];
        }        
    }

    // derive solution
    // auto juage
    if (-1!=DetIs0 && 1!=DetIs0 && 0!=DetIs0){
        printf("------[arg]DetIs0 input error!------");
        exit(EXIT_FAILURE);
    }
    if (-1 == DetIs0){
        double LambdaMin = A[0][0];
        double LambdaMax = A[0][0];
        for (int k=0;k<n;k++){
            if (fabs(A[k][k])>LambdaMax){
                LambdaMax = A[k][k];
            }
            if (fabs(A[k][k])<LambdaMin){
                LambdaMin = A[k][k];
            }
        }

        DetIs0 = (1>LambdaMin && 0.005>fabs(LambdaMin/LambdaMax)) ? 1 : 0;
    }
    // DetA != 0
    if (0 == DetIs0){
        for (int i=n-1;i>-1;i--)
        {
            for (int j=i+1;j<n;j++)
            {
                b[i] = b[i] - b[j]*A[i][j];
            }
            b[i] = b[i] / A[i][i];
        }
    }
    // DetA == 0
    else {
        if (0. != b[n-1]){
            printf("No Solution!");
            exit(EXIT_SUCCESS);
        }
        
        b[n-1] = 1.;
        for (int i=n-1;i>-1;i--)
        {
            for (int j=i+1;j<n;j++)
            {
                b[i] = b[i] - b[j]*A[i][j];
            }
            b[i] = b[i] / A[i][i];
        }
        VectNorm(n, b);        
    }

    // return X
    for (int i=0;i<n;i++){
        X[i] = b[i];
    }
}


void Gauss_Seidel_iter(int n, double **A, double *X, double *b, double epsilon)
{
    double x1[n];
    double x0[n];
    VectInit_Same(n, x1, 1);
    VectInit_Same(n, x0, 10);

    // Calculate M
    double M[n][n];
    for (int i=0;i<n;i++){   
        for (int j=0;j<n;j++){
            M[i][j] = (i==j ? 0 : (- A[i][j] / A[i][i]));
        }
    }
    // Calculate b=D^{-1}b
    for (int i=0;i<n;i++){
        b[i] = b[i] / A[i][i];
    } 

    // Iteration
    while (Distance(n, x1, x0) > epsilon){   
        // copy x1 to x0
        for (int i=0;i<n;i++){
            x0[i] = x1[i];
        }
        // update x1
        for (int i=0;i<n;i++){
            x1[i] = b[i];       
            for (int j=0;j<n;j++){
                x1[i] = x1[i] + M[i][j]*x1[j];
            }
        }
    }

    // return X
    for (int i=0;i<n;i++){
        X[i] = x1[i];
    }
}


void Eigvalue_Jacobi(int n, double **Matrix, double *Eigvalues, double epsilon)
{
    double sum_nondiag2 = 0.; // Sum(Matrix[i][j] for i != j)
    double MaxNondiag = 0.;
    int MaxNondiag_p, MaxNondiag_q;
    double s, t, cos, sin;

    // buffer
    double **MatrixBuffer;
    MatrixBuffer = MatrixGenerator(n, n, MatrixInit_0);
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            MatrixBuffer[i][j] = Matrix[i][j];
        }
    }   

    // Homework requirement, you can delete it if you think it useless.
    // Output sum_nondiag2
    // Create output file.
    char *filename = "Sum_NonDiag2.csv";
    FILE * fp = fopen(filename, "w+");
    if (NULL == fp){
        printf("fopen() failed!!!");
        exit(EXIT_FAILURE);
    }

    do
    {
        sum_nondiag2 = 0.;
        MaxNondiag = 0.;
        MaxNondiag_p = 0;
        MaxNondiag_q = 0;
        
        // 选取非对角按模最大
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                if ((i!=j) && (fabs(MaxNondiag)<fabs(MatrixBuffer[i][j]))){
                    MaxNondiag = MatrixBuffer[i][j];
                    MaxNondiag_p = i;
                    MaxNondiag_q = j;
                }
            }
        }

        // 确定旋转角theta
        s = (MatrixBuffer[MaxNondiag_q][MaxNondiag_q] - MatrixBuffer[MaxNondiag_p][MaxNondiag_p]) / 
                (2.*MatrixBuffer[MaxNondiag_p][MaxNondiag_q]);
        if (0.==s){
            t = 1.;
            cos = 1. / sqrt(2);
            sin = 1. / sqrt(2);
        }
        else {
            t = (s>0 ? (-s+sqrt(1+s*s)) : (-s-sqrt(1+s*s)));
            cos = 1. / sqrt(1+t*t);
            sin = t / sqrt(1+t*t);
        }
        

        // Update elements
        for (int i=0;i<n;i++){
            if ((i==MaxNondiag_p) || (i==MaxNondiag_q)){
                continue;
            }

            // B[i][p] = B[p][i]
            double B_ip = MatrixBuffer[MaxNondiag_p][i]*cos - MatrixBuffer[MaxNondiag_q][i]*sin;
            // B[i][q] = B[q][i]
            double B_iq = MatrixBuffer[MaxNondiag_p][i]*sin + MatrixBuffer[MaxNondiag_q][i]*cos;

            MatrixBuffer[MaxNondiag_p][i] = B_ip;
            MatrixBuffer[i][MaxNondiag_p] = B_ip;
            MatrixBuffer[MaxNondiag_q][i] = B_iq;
            MatrixBuffer[i][MaxNondiag_q] = B_iq;
        }
        double B_pp = MatrixBuffer[MaxNondiag_p][MaxNondiag_p]*cos*cos + 
            MatrixBuffer[MaxNondiag_q][MaxNondiag_q]*sin*sin - 
                MatrixBuffer[MaxNondiag_p][MaxNondiag_q]*2.*sin*cos;
        double B_qq = MatrixBuffer[MaxNondiag_p][MaxNondiag_p]*sin*sin + 
            MatrixBuffer[MaxNondiag_q][MaxNondiag_q]*cos*cos + 
                MatrixBuffer[MaxNondiag_p][MaxNondiag_q]*2.*sin*cos;
        // B[p][q] = B[q][p]
        double B_pq = MatrixBuffer[MaxNondiag_p][MaxNondiag_q]*(cos*cos - sin*sin) + 
            (MatrixBuffer[MaxNondiag_p][MaxNondiag_p]-MatrixBuffer[MaxNondiag_q][MaxNondiag_q])*sin*cos;
        MatrixBuffer[MaxNondiag_p][MaxNondiag_p] = B_pp;
        MatrixBuffer[MaxNondiag_q][MaxNondiag_q] = B_qq;
        MatrixBuffer[MaxNondiag_q][MaxNondiag_p] = B_pq;
        MatrixBuffer[MaxNondiag_p][MaxNondiag_q] = B_pq;

        // calculate Sum(Matrix[i][j] for i != j)
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                sum_nondiag2 = sum_nondiag2 + (i==j ? 0. : MatrixBuffer[i][j]*MatrixBuffer[i][j]);
            }
        }

        // Output result
        fprintf(fp, "%f\n", sum_nondiag2);

    } while (sum_nondiag2 > epsilon);
    printf((0 == fclose(fp)) ? "File closed successfully.\n" : "Failed to close file."); 


    // return results
    printf("EigenValues: \n");
    for (int i=0;i<n;i++){
        Eigvalues[i] = MatrixBuffer[i][i];
        // homework requirement.
        printf("%f  ", Eigvalues[i]);
    }
    printf("\n");

    FreeMatrixMemory(MatrixBuffer, n);
}


void Eigvector(int n, double **Matrix, double *Eigvectors, double EigValue, 
                void(*solver)(int, double**, double*, double*, int))
{
    double b[n];
    VectInit_Same(n, b, 0);

    // buffer
    double **A;
    A = MatrixGenerator(n, n, MatrixInit_0);
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            A[i][j] = (i==j ? Matrix[i][j]-EigValue : Matrix[i][j]);
        }
    }   

    solver(n, A, Eigvectors, b, 1);

    FreeMatrixMemory(A, n);
}


void Eigvector_iter(int n, double **Matrix, double *Eigvectors, double EigValue, double epsilon,
                void(*solver)(int, double**, double*, double*, double))
{
    double b[n];
    VectInit_Same(n, b, 0);

    // buffer
    double **A;
    A = MatrixGenerator(n, n, MatrixInit_0);
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            A[i][j] = (i==j ? Matrix[i][j]-EigValue : Matrix[i][j]);
        }
    }   

    solver(n, A, Eigvectors, b, epsilon);

    FreeMatrixMemory(A, n);
}


void SVD(int n, int m, double **A, double **U, double **Sigma, 
            double **VT, double epsilon)
{
    double **AT = MatrixGenerator(m, n, MatrixInit_0);
    double **A_AT = MatrixGenerator(n, n, MatrixInit_0);
    double **AT_A = MatrixGenerator(m, m, MatrixInit_0);
    double A_AT_eigvalues[n];
    double AT_A_eigvalues[m];

    MatrixTranspose(n, m, A, AT);
    MatrixMutiply(n, m, A, AT, A_AT);
    MatrixMutiply(m, n, AT, A, AT_A);


    // Homework requirement. You can delete it if you think it useless.
    printf("A*A^T : \n");
    for (int i=0;i<n;i++){
        printf("{ ");
        for (int j=0;j<n;j++){
            printf("%f  ", A_AT[i][j]);            
        }
        printf("\b}\n");
    }


    // A*AT => U
    Eigvalue_Jacobi(n, A_AT, A_AT_eigvalues, epsilon);    
    // sort eigenvalues >
    for (int i=0;i<n;i++){
        double term = 0.;
        for (int j=i+1;j<n;j++){
            if (A_AT_eigvalues[j] > A_AT_eigvalues[i]){
                term = A_AT_eigvalues[j];
                A_AT_eigvalues[j] = A_AT_eigvalues[i];
                A_AT_eigvalues[i] = term;
            }
        }
    }
    // fill U
    for (int i=0;i<n;i++){
        double EigVect[n]; 
        Eigvector(n, A_AT, EigVect, A_AT_eigvalues[i], GaussElim_PartialPivoting);
        for (int j=0;j<n;j++){
            U[j][i] = EigVect[j];
        }
    }


    // AT*A => V
    Eigvalue_Jacobi(m, AT_A, AT_A_eigvalues, epsilon);
    // sort eigenvalues
    for (int i=0;i<m;i++){
        double term = 0.;
        for (int j=i+1;j<m;j++){
            if (AT_A_eigvalues[j] > AT_A_eigvalues[i]){
                term = AT_A_eigvalues[j];
                AT_A_eigvalues[j] = AT_A_eigvalues[i];
                AT_A_eigvalues[i] = term;
            }
        }
    }
    // fill VT
    for (int i=0;i<m;i++){
        double EigVect[m]; 
        Eigvector(m, AT_A, EigVect, AT_A_eigvalues[i], GaussElim_PartialPivoting);
        for (int j=0;j<m;j++){
            VT[i][j] = EigVect[j];
        }
    }

    // fill Sigma
    MatrixInit_0(n, m, Sigma);
    double Min_n_m = (m<n ? m : n);
    for (int i=0;i<Min_n_m;i++){
        Sigma[i][i] = sqrt(A_AT_eigvalues[i]);
    }


    FreeMatrixMemory(AT, m);
    FreeMatrixMemory(A_AT, n);
    FreeMatrixMemory(AT_A, m);
}


void PCA(int n, int m, double **M, double epsilon, int FinalDim, double **FinalM, int NeedBUffer)
{
    // buffer
    double **A;
    if (1==NeedBUffer){
        A = MatrixGenerator(n, m, MatrixInit_0);
        for (int i=0;i<n;i++){
            for (int j=0;j<m;j++){
                A[i][j] = M[i][j];
            }
        }       
    }
    else {
        A = M;
    }

    // uncerternize
    double meanVect[n];
    VectInit_Same(n, meanVect, 0);
    for (int j=0;j<m;j++){
        for (int i=0;i<n;i++){
            meanVect[i] = meanVect[i] + A[i][j];
        }
    }
    for (int i=0;i<n;i++){
        meanVect[i] = meanVect[i] / m;
    }
    for (int j=0;j<m;j++){
        for (int i=0;i<n;i++){
            A[i][j] = A[i][j] - meanVect[i];
        }        
    }


    double **AT = MatrixGenerator(m, n, MatrixInit_0);
    double **A_AT = MatrixGenerator(n, n, MatrixInit_0);
    double A_AT_eigvalues[n];
 
    MatrixTranspose(n, m, A, AT);
    MatrixMutiply(n, m, A, AT, A_AT);


    // Homework requirement. You can delete it if you think it useless.
    printf("X*X^T/m : \n");
    for (int i=0;i<n;i++){
        printf("[ ");
        for (int j=0;j<n;j++){
            printf("%f  ", A_AT[i][j]/m);            
        }
        printf("\b]\n");
    }


    Eigvalue_Jacobi(n, A_AT, A_AT_eigvalues, epsilon);    
    // sort eigenvalues >
    for (int i=0;i<n;i++){
        double term = 0.;
        for (int j=i+1;j<n;j++){
            if (A_AT_eigvalues[j] > A_AT_eigvalues[i]){
                term = A_AT_eigvalues[j];
                A_AT_eigvalues[j] = A_AT_eigvalues[i];
                A_AT_eigvalues[i] = term;
            }
        }
    }

    // Calculate FinalM
    for (int i=0;i<FinalDim;i++){
        double ProjectionBase[n]; 
        Eigvector(n, A_AT, ProjectionBase, A_AT_eigvalues[i], GaussElim_PartialPivoting);
        for (int j=0;j<m;j++){
            FinalM[i][j] = 0.;
            for (int k=0;k<n;k++){
                FinalM[i][j] = FinalM[i][j] + A[k][j]*ProjectionBase[k];
            }
        }
    }


    FreeMatrixMemory(A, n);
    FreeMatrixMemory(AT, m);
    FreeMatrixMemory(A_AT, n);
}


