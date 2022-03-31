#include <stdio.h>
#define SPACE -12

double dAbs(double num);
void initRangeMat(int n, int num[n]);

void printNMat(int rows, int cols, double mat[rows][cols], int width);
void printOneMat(int rows, int cols, double *mat, int width);
void printRangeMat(int rows, int cols, int *mat, int width);

void rangeRowsNMat(int n, double mat[n][n], int row1, int row2);
void rangeRolsOneMat(int cols, double mat[cols], int col1, int col2);
void rangeRolsRangeMat(int cols, int mat[cols], int col1, int col2);

void gauss(int n, double mat[n][n], double b[n], double x[n], int x_range[n], int main_x);

void lu(int n, double mat[n][n]);
void solveLMat(int n, double l[n][n], double y[n], double b[n]);
void solveUMat(int n, double u[n][n], double x[n], double y[n]);
void solveLU(int n, double mat[n][n], double b[n], double result[n]);

int main(void)
{
    double num[4][4] = {{1, 2, 3, -4}, {-3, -4, -12, 13}, {2, 10, 0, -3}, {4, 14, 9, -13}};
    double b[4] = {-2, 5, 10, 7};
    double x[4];
    double y[4];
    int x_range[4];

    solveLU(4, num, b, x);
    printOneMat(1, 4, x, SPACE);

    return 0;
}

double dAbs(double num)
{
    if (num >= 0)
        return num;
    else
        return -num;
}

void initRangeMat(int n, int num[n])
{
    for (int i = 0; i < n; i++)
        num[i] = i + 1;
}

void printNMat(int rows, int cols, double mat[rows][cols], int width)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
        {
            printf("% *g", width, mat[i][j]);
            if (j == cols - 1)
                printf("\n");
        }
    printf("\n");
}

void printOneMat(int rows, int cols, double *mat, int width)
{
    if (rows == 1)
    {
        for (int i = 0; i < cols; i++)
            printf("% *g", width, *mat++);
    }
    else if (cols == 1)
        for (int i = 0; i < rows; i++)
            printf("% *g\n", width, *mat++);
    printf("\n");
}

void rangeRowsNMat(int n, double mat[n][n], int row1, int row2)
{
    double temp[n];
    for (int i = 0; i < n; i++)
    {
        temp[i] = mat[row1][i];
        mat[row1][i] = mat[row2][i];
        mat[row2][i] = temp[i];
    }
}

void rangeRolsOneMat(int cols, double mat[cols], int col1, int col2)
{
    double temp;
    for (int i = 0; i < cols; i++)
    {
        temp = mat[col1];
        mat[col1] = mat[col2];
        mat[col2] = temp;
    }
}

void rangeRolsRangeMat(int cols, int mat[cols], int col1, int col2)
{
    int temp;
    for (int i = 0; i < cols; i++)
    {
        temp = mat[col1];
        mat[col1] = mat[col2];
        mat[col2] = temp;
    }
}

void printRangeMat(int rows, int cols, int *mat, int width)
{
    if (rows == 1)
    {
        for (int i = 0; i < cols; i++)
            printf("x%*d", width, *mat++);
    }
    else if (cols == 1)
        for (int i = 0; i < rows; i++)
            printf("x%*d", width, *mat++);
    printf("\n");
}

void gauss(int n, double mat[n][n], double b[n], double x[n], int x_range[n], int main_x)
{
    // n是n阶矩阵,s是对角线元素,r是行,c是列,k是倍数,x_range是换行之后的排列,main_x是高斯列主元消去法
    int rows, cols;
    rows = cols = n;
    double k, max;
    for (int s = 0; s < rows - 1; s++)
    {
        max = dAbs(mat[s][s]);
        for (int r = s + 1; r < rows; r++)
        {
            if (main_x)
            {
                if (dAbs(mat[r][s]) > max)
                {
                    max = mat[r][s];
                    rangeRolsOneMat(n, b, r, s);
                    rangeRolsRangeMat(n, x_range, r, s);
                    rangeRowsNMat(n, mat, r, s);
                }
            }

            k = -mat[r][s] / mat[s][s];
            b[r] += b[s] * k;

            for (int c = 0; c < rows; c++)
                mat[r][c] += mat[s][c] * k;
        }
    }

    x[rows - 1] = b[rows - 1] / mat[rows - 1][rows - 1];
    for (int s = rows - 2; s >= 0; s--)
    {
        for (int cols = rows - 1; cols > s; cols--)
            b[s] -= mat[s][cols] * x[cols];
        x[s] = b[s] / mat[s][s];
    }
}

void lu(int n, double mat[n][n])
{
    for (int r = 1; r < n; r++)
        mat[r][0] /= mat[0][0];
    for (int k = 1; k < n; k++) // k是对角线元素所在的行(列)
    {
        for (int r = k; r < n; r++)
            for (int i = 0; i < k; i++)             // i是从第0行(列)迭代到对角线元素
                mat[k][r] -= mat[k][i] * mat[i][r]; //向右-=

        for (int c = k + 1; c < n; c++)
        {
            for (int j = 0; j < k; j++)             // j是从第0行(列)迭代到对角线元素
                mat[c][k] -= mat[c][j] * mat[j][k]; //向下-/
            mat[c][k] /= mat[k][k];
        }
    }
}

void solveLMat(int n, double l[n][n], double y[n], double b[n])
{
    double temp;
    y[0] = b[0];
    for (int s = 1; s < n; s++)
    {
        temp = 0;
        for (int i = 0; i < s; i++)
            temp += l[s][i] * y[i];
        y[s] = b[s] - temp;
    }
}

void solveUMat(int n, double u[n][n], double x[n], double y[n])
{
    double temp;
    x[n - 1] = y[n - 1] / u[n - 1][n - 1];
    for (int s = n - 2; s >= 0; s--)
    {
        temp = 0;
        for (int i = n - 1; i > s; i--)
            temp += u[s][i] * x[i];
        x[s] = (y[s] - temp) / u[s][s];
    }
}

void solveLU(int n, double mat[n][n], double b[n], double result[n])
{
    double y[n];
    lu(n, mat);
    solveLMat(n, mat, y, b);
    solveUMat(n, mat, result, y);
}
