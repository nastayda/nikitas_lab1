//
//  main.c
//  Lab_1
//
//  Created by Анастасия Данилкина on 21.11.17.
//  Copyright © 2017 Анастасия Данилкина. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

struct result_t{
    std::size_t countOfIterations;
    double error;
}hola;

double Sqr (const double x ){
    return x * x;
}

void TridiagonalMatrixAlgorithm ( const int size,
                                 const std::vector< std::vector< double >> & Diagonals,
                                 const std::vector< double > & vectorB,
                                 std::vector< double > & vectorResult ){
    double *alpha = new double [size];
    double *beta = new double [size];
    double *gamma = new double [size];
    gamma[0] = Diagonals[0][1];
    beta[0] = vectorB[0] / gamma[0];
    alpha[0] = - Diagonals[0][2] / gamma[0];
    for(int i=1;i<size;i++){
        gamma[i] = Diagonals[i][1] + Diagonals[i][0] * alpha[i - 1];
        beta[i] = ( vectorB[i] - Diagonals[i][0] * beta[i - 1] ) / gamma[i];
        alpha[i] = - Diagonals[i][2] / gamma[i];
    }
    vectorResult[size - 1] = beta [size - 1];
    for(int i=size-2;i>=0;i--){
        vectorResult[i] = alpha[i] * vectorResult[i + 1] + beta[i];
    }
    delete []alpha;
    delete []beta;
    delete []gamma;
}

void InitialBoundaryValueProblem (
                                  std::vector< std::vector< double >> & x,
                                  const std::vector< std::vector< double >> & f,
                                  const std::vector< double > & p,
                                  const double a,
                                  const double l,
                                  const double v, const double T, const int lCount, const int tCount ){
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    const double h=l/(lCount-1.0);
    const double tau=T/(tCount-1.0);
    const double q=0.5*h*h/(a*a);
    A.resize( lCount );
    b.resize( lCount );
    //Генерируем матрицу системы
    for(int i=0;i<lCount;i++){ A[i].resize( 3 );
        if ( i == 0 ) {
            A[0][1] = tau + q;
            A[0][2] = -tau;
        }else {
            if ( i == lCount - 1 ) {
                A[lCount - 1][1] = tau + tau * h * v + q;
                A[lCount - 1][0] = -tau;
            } else {
                A[i][1] = h * h + 2.0 * a * a * tau;
                A[i][0] = -a * a * tau;
                A[i][2] = -a * a * tau;
            }
        }
    }
    for(int i=1;i<tCount;i++){
        //Заполняем массив правых частей
        for ( int j = 0; j < lCount; j++ ) {
            if ( j == 0 ) {
                b[j] = q * ( tau * f[i - 1][j] + x[i - 1][j] );
            } else {
                if ( j == lCount - 1 ) {
                    b[j] = tau * h * v * p[i - 1] + q * ( tau * f[i - 1][j] + x[i - 1][j] );
                } else {
                    b[j] = h * h * ( x[i - 1][j] + tau * f[i - 1][j] ); }
            }
        }
        //Решаем систему
        TridiagonalMatrixAlgorithm( lCount, A, b, x[i] ); }
}


//NE NUZHEN
double NormL2 (
               const double h,
               const double tau,
               const std::vector< std::vector< double >> f ){
    double result = 0;
    const std::size_t N = f.size(), M = f[0].size();
    //Внутренность по i
    for(auto i=1;i<N-1;i++){
        //Внутренность по j
        for ( auto j = 1; j < M - 1; j++ ) {
            result += Sqr( f[i][j] );
        }
        //Края по j
        result += 0.5 * ( Sqr( f[i][0] ) + Sqr( f[i][M - 1] ) );
    }
     //Края по i
     //i==0
     for(auto j=1;j<M-1;j++){
     //Внутренность по j
     result += 0.5 * Sqr( f[0][j] );
     }
     //Края по j
     result += 0.25 * ( Sqr( f[0][0] ) + Sqr( f[0][M - 1] ) );
     //Края по i
     //i==N-1
     for(auto j=1;j<M-1;j++){
     //Внутренность по j
     result += 0.5 * Sqr( f[N - 1][j] );
    }
    //Края по j
    result+=0.25*(Sqr(f[N-1][0])+Sqr(f[N-1][M-1]));
    //Итог
    result = std::sqrt( h * tau * result );
    return result;
}

//DRUGOI!!!
void Projection (
                 std::vector< std::vector< double >> & f,
                 const std::vector< std::vector< double >> & ph, const double R,
                 const double h,
                 const double tau ){
    const double _1_norm = 1.0 / NormL2( h, tau, ph );
    std::size_t N = ph.size(), M = ph[0].size();
    for(auto i=0;i<N-1;i++){
        for ( auto j = 0; j < M; j++ ) {
            f[i][j] = -R * _1_norm * ph[i + 1][j]; }
    }
}

void AdditionalInitialBoundaryValueProblem ( std::vector<std::vector< double >> & ph,
                                            const double a,
                                            const double l,
                                            const double v,
                                            const double T,
                                            const int lCount,
                                            const int tCount ){
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    const double h=l/(lCount-1.0);
    const double tau=T/(tCount-1.0);
    const double q=0.5*h*h/(a*a);
    A.resize( lCount );
    b.resize( lCount );
    //Генерируем матрицу системы
    for(int i=0;i<lCount;i++){
        A[i].resize( 3 );
        if ( i == 0 ) {
            A[0][1] = tau + q;
            A[0][2] = -tau; }else {
                if ( i == lCount - 1 ) {
                    A[lCount - 1][1] = tau + tau * h * v + q;
                    A[lCount - 1][0] = -tau;
                } else {
                    A[i][1] = ( h * h + 2.0 * a * a * tau );
                    A[i][0] = -a * a * tau;
                    A[i][2] = -a * a * tau;
                } }
    }
    //Итерации по времени
    for(int i=tCount-2;i>=0;i--){
        //Заполняем массив правых частей
        for ( int j = 0; j < lCount; j++ ) {
            if ( j == 0 ) {
                b[j] = q * ( ph[i + 1][j] );
            } else {
                if ( j == lCount - 1 ) {
                    b[j] = q * ( ph[i + 1][j] );
                }
                else {
                    b[j] = h * h * ( ph[i + 1][j] );
                }
            }
        }
        //Решаем систему
        TridiagonalMatrixAlgorithm( lCount, A, b, ph[i] ); }
}

double CalculateAlpha1 (
                        const double tau,
                        const std::vector< std::vector< double >> & ps, const std::vector< std::vector< double >> & f, const std::vector< std::vector< double >> & f_,
                        const std::vector< std::vector< double >> & x, const std::vector< std::vector< double >> & x_ ){
    double alpha = 1.0;
    std::size_t N = f.size(), M = f[0].size();
    double result = 0;
    double buf = 0;
    //Вычисляем числитель
    //Внутренность по i
    for(auto i=1;i<N-1;i++){
        //Внутренность по j
        for ( auto j = 1; j < M - 1; j++ ) {
            result += ps[i + 1][j] * ( f_[i][j] - f[i][j] );
        }
        //Края по j
        result += 0.5 * ( ps[i + 1][0] * ( f_[i][0] - f[i][0] ) + ps[i + 1][M - 1] * ( f_[i][M - 1] - f[i][M - 1] ) );
    }
    //Края по i
    //i==0
    for(auto j=1;j<M-1;j++){
        //Внутренность по j
        result += 0.5 * ps[1][j] * ( f_[0][j] - f[0][j] );
    }
    //Края по j
    result+=0.25*(ps[1][0]*(f_[0][0]-f[0][0])+ps[1][M-1]*(f_[0][M-1]-f[0][M-1]));
    //Края по i
    //i==N-1
    for(auto j=1;j<M-1;j++){
        //Внутренность по j
        result += 0.5 * ps[N][j] * ( f_[N - 1][j] - f[N - 1][j] );
    }
    //Края по j
    result+=0.25*(ps[N][0]*(f_[N-1][0]-f[N-1][0])+ps[N][M-1]*(f_[N-1][M-1]-f[N- 1][M - 1] ) );
    //Итог
    result = tau * result;
    //Вычисляем знаменатель
    for(auto i=1;i<M-1;i++){
        buf += 2 * Sqr( x_[N][i] - x[N][i] );
    }
    buf += ( Sqr( x_[N][0] - x[N][0] ) + Sqr( x_[N][M - 1] - x[N][M - 1] ) );
    //Вычисляем и возвращаем нужное
    alpha = - result / buf;
    if(alpha<1.0){
        if ( alpha < 0 ) {
            return 0;
        }
        return alpha;
    } else {
        return 1.0;
    }
}


double ConditionalGradientMethodStep (
                                      std::vector< std::vector< double >> & f,
                                      const std::vector< std::vector< double >> & f_,
                                      const double alpha,
                                      const double R,
                                      const double h,
                                      const double tau ){
    std::size_t N = f.size(), M = f[0].size();
    double result = 0, buf = 0;
    //Внутренность по i
    for(auto i=1;i<N-1;i++){
        //Внутренность по j
        for ( auto j = 1; j < M - 1; j++ ) {
            buf = f[i][j];
            f[i][j] = ( 1.0 - alpha ) * f[i][j] + alpha * f_[i][j];
            buf = buf - f[i][j];
            result += Sqr( buf );
        }
        //Край j = 0
        buf = f[i][0];
        f[i][0] = ( 1.0 - alpha ) * f[i][0] + alpha * f_[i][0];
        buf = buf - f[i][0];
        result += 0.5 * Sqr( buf );
        //Край j = M - 1
        buf = f[i][M - 1];
        f[i][M - 1] = ( 1.0 - alpha ) * f[i][M - 1] + alpha * f_[i][M - 1];
        buf= buf-f[i][M-1];
        result += 0.5 * Sqr( buf );
    }
    //Края по i
    //i==0
    //Внутренность по j
    for(auto j=1;j<M-1;j++){
        buf = f[0][j];
        f[0][j] = ( 1.0 - alpha ) * f[0][j] + alpha * f_[0][j];
        buf = buf - f[0][j];
        result += 0.5 * Sqr( buf );
    }
    //Крайj=0
    buf = f[0][0];
    f[0][0] = ( 1.0 - alpha ) * f[0][0] + alpha * f_[0][0];
    buf = buf - f[0][0];
    result += 0.25 * Sqr( buf );
    //Крайj=M-1
    buf = f[0][M - 1];
    f[0][M-1]=(1.0-alpha)*f[0][M-1]+alpha*f_[0][M-1];
    buf=buf-f[0][M-1];
    result += 0.25 * Sqr( buf );
    //i==N-1
    //Внутренность по j
    for(auto j=1;j<M-1;j++){
        buf = f[N - 1][j];
        f[N - 1][j] = f[N - 1][j] + alpha * ( f_[N - 1][j] - f[N - 1][j] );
        buf= buf-f[N-1][j];
        result += Sqr( buf );
    }
    //Крайj=0
    buf = f[N - 1][0];
    f[N-1][0]=(1.0-alpha)*f[N-1][0]+alpha*f_[N-1][0];
    buf=buf-f[N-1][0];
    result += 0.5 * Sqr( buf );
    //Крайj=M-1
    buf=f[N-1][M-1];
    f[N-1][M-1]=(1.0-alpha)*f[N-1][M-1]+alpha*f_[N-1][M-1];
    buf= buf-f[N-1][M-1];
    result += 0.5 * Sqr( buf );
    //Итог
    result = std::sqrt( h * tau * result );
    return result;
}


result_t ConditionalGradientMethod1 ( const int lCount,
                                     const int tCount,
                                     const bool save,
                                     const double a,
                                     const double l,
                                     const double v,
                                     const double T,
                                     const double R,
                                     const double eps,
                                     const std::vector< double > & p,
                                     const std::vector< double > & y,
                                     std::vector< std::vector< double >> & x,
                                     std::vector< std::vector< double >> & f ){
    double d=1.0,tau=T/(tCount-1.0),h = l / ( lCount - 1.0 ), alpha = 0;
    result_t result;
    std::vector<std::vector<double>> ps;
    std::vector<std::vector<double>> f_;
    std::vector<std::vector<double>> x_ = x;
    std::ofstream ofs( "met1.txt" );
    std::ofstream ofsF( "met1_f.txt" );
    if(save){
        for ( auto i = 0; i < lCount; i++ ) {
            ofs << i * h << " ";
        }
        ofs << std::endl;
        
        for ( auto i = 0; i < lCount; i++ ) {
            ofs << y[i] << " ";
        }
        ofs << std::endl;
    }
    ps.resize( tCount );
    f_.resize( tCount - 1 );
    for(int i=0;i<tCount-1;i++){
        ps[i].resize( lCount );
        f_[i].resize( lCount );
    }
    ps[tCount - 1].resize( lCount );
    for ( result.countOfIterations = 1; d > eps; result.countOfIterations++) {
        //Вычисляем функцию psi
        for ( int i = 0; i < lCount; i++ ) {
            ps[tCount - 1][i] = 2.0 * ( x[tCount - 1][i] - y[i] );
            if(i==0){
                result.error = 0.5 * Sqr( x[tCount- 1][i] - y[i] );
            }else{
                if ( i == lCount - 1 ) {
                    result.error += 0.5 * Sqr( x[tCount - 1][i] - y[i] );
                } else {
                    result.error += Sqr( x[tCount - 1][i] - y[i] );
                }
            }
        }
        result.error = std::sqrt( h * result.error );
        AdditionalInitialBoundaryValueProblem( ps, a, l, v, T, lCount, tCount );
        //Находим вспомогательное приближение
        Projection( f_, ps, R, h, tau );
        //Вычисляем alpha
        InitialBoundaryValueProblem( x_, f_, p, a, l, v, T, lCount, tCount );
        alpha = CalculateAlpha1( tau, ps, f, f_, x, x_ );
        //Вычисляем следующее значение управления и определяем норму разности с предыдущим значением
        d = ConditionalGradientMethodStep( f, f_, alpha, R, h, tau );
        //Вычисляем следующее приближение к решению краевой задачи
        InitialBoundaryValueProblem( x, f, p, a, l, v, T, lCount, tCount );
        //Сохраняем результат
        if ( save ) {
            for ( auto i = 0; i < f.size(); i++ ) {
                for ( auto j = 0; j < f[0].size(); j++ ) {
                    ofsF << std::fixed << std::setprecision( 8 ) << f[i][j] << " ";
                }
                ofsF << std::endl; }
            switch ( result.countOfIterations) {
                case 1:
                    for ( auto i = 0; i < lCount; i++ ) {
                        ofs << x[tCount - 1][i] << " ";
                    }
                    ofs << std::endl; break;
                case 2:
                    for ( auto i = 0; i < lCount; i++ ) {
                        ofs << x[tCount - 1][i] << " ";
                    }
                    ofs << std::endl;
                    break; case 3:
                    for ( auto i = 0; i < lCount; i++ ) {
                        ofs << x[tCount - 1][i] << " ";
                    }
                    ofs << std::endl; break;
                case 4:
                    for ( auto i = 0; i < lCount; i++ ) {
                        ofs << x[tCount - 1][i] << " ";
                    }
                    ofs << std::endl;
                    break; case 6:
                    for ( auto i = 0; i < lCount; i++ ) {
                        ofs << x[tCount - 1][i] << " ";
                    }
                    ofs << std::endl;
                    break; case 8:
                    for ( auto i = 0; i < lCount; i++ ) {
                        ofs << x[tCount - 1][i] << " ";
                    }
                    ofs << std::endl; break;
                case 10:
                    for ( auto i = 0; i < lCount; i++ ) {
                        ofs << x[tCount - 1][i] << " ";
                    }
                    ofs << std::endl;
                    break; case 100:
                    for ( auto i = 0; i < lCount; i++ ) {
                        ofs << x[tCount - 1][i] << " ";
                    }
                    ofs << std::endl; break;
                case 500:
                    for ( auto i = 0; i < lCount; i++ ) {
                        ofs << x[tCount - 1][i] << " ";
                    }
                    ofs << std::endl;
                    break;
            }
        }
    }
    return result;
}

int main(int argc, const char * argv[]) {
    
    std::vector<std::vector<double>> x1, f1, x2, f2;
    std::vector<double> p;
    std::vector<double> y;
    const unsigned int precision = 10;
    const int lCount = 100.0;
    const int tCount = 200.0;
    const double R = 30;
    const double eps = std::pow( 10.0, - (double) precision );
    const double a = 5.0;
    const double l = 2.0;
    const double T = 3.0;
    const double v = 20.0;
    const double h=l/(lCount-1.0);
    const double tau=T/(tCount-1.0);
    std::ofstream ofs( "res.txt" );
    y.resize( lCount );
    p.resize( tCount - 1 );
    x1.resize( tCount );
    f1.resize( tCount - 1 );
    x2.resize( tCount );
    f2.resize( tCount - 1 );
    for(int i=0;i<tCount;i++){
        if ( i == 0 ) {
            x1[i].resize( lCount );
            x2[i].resize( lCount );
            for ( int j = 0; j < lCount; j++ ) {
                x1[i][j] = 4.0;
                x2[i][j] = 4.0;
            }
        }else {
            p[i - 1] = sin( i * tau ) * ( cos( l ) - sin( l ) / v ) + 4.0; x1[i].resize( lCount );
            f1[i - 1].resize( lCount );
            x2[i].resize( lCount );
            f2[i - 1].resize( lCount );
            for ( int j = 0; j < lCount; j++ ) {
                if ( i == tCount - 1 ) {
                    y[j] = sin( T ) * cos( 10 * j * h ) + 4;
                }
            }
        }
    }
    result_t res = ConditionalGradientMethod1( lCount, tCount, true, a, l , v, T, R, eps, p, y, x1, f1 );
    std::cout << "Точность: " << std::fixed << std::setprecision( precision ) << eps << std::endl;
    std::cout << "Метод 1" << std::endl;
    std::cout << " " << "Количество итераций: " << res.countOfIterations << std::endl;
    std::cout << " " << "Погрешность решения: " << std::fixed << std::setprecision( precision ) << res.error << std::endl;

    /*res = ConditionalGradientMethod2( lCount, tCount, true, a, l, v, T, R, eps, p, y, x2, f2 );
     std::cout << "Метод 2" << std::endl;
     std::cout << " " << "Количество итераций: " << res.countOfIterations << std::endl;
     std::cout << " " << "Погрешность решения: " << std::fixed << std::setprecision( precision ) << res.error << std::endl;*/
    
    for(auto i=0;i<lCount;i++){
        ofs << i * h << " ";
    }
    ofs << std::endl;
    
    for(auto i=0;i<lCount;i++){
        ofs << y[i] << " ";
    }
    ofs << std::endl;
    
    for(auto i=0;i<lCount;i++){
        ofs << x1[tCount - 1][i] << " ";
    }
    ofs << std::endl;
    
    for(auto i=0;i<lCount;i++){
        ofs << x2[tCount - 1][i] << " ";
    }
    ofs << std::endl;
    
    return 0;
}
