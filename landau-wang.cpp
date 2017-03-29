/*
 * landau-wang.cpp
 *
 *  Created on: 24 мая 2014 г.
 *      Author: nlare
 */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <omp.h>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <boost/filesystem.hpp>
#include "mersenne.cpp"

#define DEBUG
#define energy(b) (2*(b)-2.0*(L*L))
//#define magnet(b) (2*(b)-2.0*(L*L))

int **spin;     // Массив спинов
int L;          // Размер решетки

int Lx, Ly, Lz; // Решетка для Гейзенберга 

int neighbour_spins(int,int);
double NeighboursHeisenberg(double***,int,int,int);

int LandauWangSerialIsing();
int LandauWangSerialHeisenberg();
int DataProcessing();

int main(int argc, char *argv[])  {

    // LandauWangSerialHeisenberg();
    LandauWangSerialIsing();
    // DataProcessing();

    return 0;

}

int LandauWangSerialIsing() {

    int *hist;              // Массив гистограммы энергий, т.е количества посещений данного энергетического состояния
    double *g;              /* Массив плотности состояний
                             * Изначально каждый элемент массива принимается равным единице
                             */

    double *m, *m_sum;      /* Массив для намагниченности
                             * Изначально каждый элемент массива принимается равным нулю
                             */
    
    double *m2, *m2_sum;    /* Массив для квадрата намагниченности
                             * Изначально каждый элемент массива принимается равным нулю
                             */

    int b, b_new, top_b;    // Текущий энергетический уровень и последующий, top_b - предельный энергетический уровень
    double f, f_min, ln_f;  /* "f" - начальный множитель для энергетических уровней, 
                             * "f_min" - минимальное значение множителя
                             * на каждый принятый шаг f_m = (f_m-1)^1/2
                             */
    int skip;                // Количество шагов с неизменной энергией
    int mcs;
    int min_steps;
    double prob;             // Вероятность изменения энергетического уровня
    double flat_threshold;   // Порог "плоскости" гистограммы

    double time_b, time_e;

    double T;
    double seed;

    int it_count;
    
    std::stringstream ss;
    std::ofstream test_g_f, graph_g_f, graph_sh_g;

    ss.str("");
    ss << "test_g";
    boost::filesystem::create_directories(ss.str().c_str());

    srand(time(NULL));
    seed = 1 + rand() % 10000;

    CRandomMersenne Mersenne(seed);

    L = 48;
    f = 2.7182818284;  // В работе Ландау-Ванга было указано значение "f" равное экспоненте 
    f_min = 1.0000001;  // Данная переменная должна быть около единицы
    min_steps = 10000;
    skip = 10000;
    flat_threshold = 0.8;
    it_count = 0;

    // L *= 0.75; 

    // top_b=L*L*0.75;

    spin = new int* [L];
    for(int i = 0; i < L; i++)    {
        spin[i] = new int [L];
    };
    hist = new int [4*L*L];
    g = new double [4*L*L];

    m = new double [4*L*L];
    m_sum = new double [4*L*L];

    m2 = new double [4*L*L];
    m2_sum = new double [4*L*L];
    
    for(int i = 0; i < L; i++)    {
        for(int j = 0; j < L; j++)    {
            spin[i][j] = 1; 
        }
    }

    for(int i = 0; i < 2*L*L; i++)    {
        g[i] = 1.0;
    }

    for(int i = 0; i < 2*L*L; i++)    {
        m[i] = 0.0;
        m_sum[i] = 0.0;
    }

    for(int i = 0; i < 2*L*L; i++)    {
        m2[i] = 0.0;
        m2_sum[i] = 0.0;
    }

    for(int i = 0; i < 2*L*L; i++)    {
        hist[i] = 1.0;
    }

    if(L%2!=0) top_b=2*(L*(L-1))+1;
    else top_b=2*(L*L)+1;

    // top_b *= 0.75; 

    b = 0;

    #ifdef DEBUG
    std::cout << std::setprecision(10) \
              << "Energy level range: top_b = " << top_b << std::endl \
              << "Precision: f_min = " << f_min << std::endl; 
    #endif

    time_b = omp_get_wtime();

    // Считаем пока "f" не примет минимальное значение
    while(f > f_min)    {

        int count, n;
        int c = skip+1;

        ln_f = log(f);  // Для вычислений будем использовать логарифм множителя "f" 
        count = 1;  // Если гистограмма не стала плоской, останется равным единице и цикл продолжится.
        mcs = 0;
        n = 0;

        for (int i = 0; i < 4*L*L; i++)   {
            hist[i] = 0;    
        }

        // std::cout << "chck = " << b << std::endl;

        do {

        for(int i = 0; i < L; i++)   {
            for(int j = 0; j < L; j++)  {
                m1:
                int ci = Mersenne.IRandomX(0, L-1);
                int cj = Mersenne.IRandomX(0, L-1);

                if(spin[ci][cj] == 0) goto m1;  // Исключаем немагнитные спины

                b_new = b + spin[ci][cj]*neighbour_spins(ci,cj);   // Считаем новое состояние

                if(b_new < top_b)   {   // Если энергия в допустимых пределах
                    
                    prob = exp(g[b]-g[b_new]);  // Вероятность принятия нового энергетического состояния

                    if((double)(Mersenne.IRandomX(1,10000)/10000.0) < prob)    {

                        b = b_new;
                        spin[ci][cj] *= -1;

                    }
                    
                    g[b] += ln_f;         // Увеличиваем текущий энергетический уровень на логарифм "f"

                    // for(int i = 0; i < top_b; i++)  {

                    m[b] += spin[ci][cj];                 // Увеличиваем значение намагниченности для данного энергетического уровня

                    m2[b] += spin[ci][cj]*spin[ci][cj];   // Квадрат намагниченности для расчета восприимчивости

                    // }

                    hist[b] += 1;
                    n++;
                }


            }
        }

        mcs++;
        c++;
        
        if((mcs >= min_steps) && (c >= skip))    {
            
            c = 0;
        
            #ifdef DEBUG
            // std::cout << "In count area. \n";
            // std::cout << "mcs = " << mcs << std::endl;
            #endif

            int h_count = 0;
            double h_delt;
            double h_sum = 0.0;

            for(int i = 0; i < top_b; i++) {

                if(hist[i] != 0)    {
                    h_sum += (double)hist[i]/min_steps;
                    h_count++;
                }

            }
            // for(int i = 0; i < top_b; i++)
            // std::cout << "g[" << i << "]=" << g[i] << std::endl;

            count = 0;  // Выходим из цикла по окончанию

            for (int i = 0; i < top_b; i++)    {
                if(hist[i] != 0)    {
                    h_delt = (double)hist[i]/(h_sum/h_count)/min_steps;
                    if((h_delt <= flat_threshold) || (h_delt >= 1.0+flat_threshold)) count = 1;
                }
            } 

            std::cout << "ln_f=" << ln_f << ", flat_threshold = " << flat_threshold << ", h_delt = " << h_delt << ", count = " << count << std::endl;

        }   else count = 1;

        }   while(count);

        for(int i = 1; i <= top_b; i++) {

            g[i] -= g[0]; 

            // m_sum[i] = 0.0;
            m_sum[i] += m[i];
            m[i] = 0.0;
            
            // m2_sum[i] = 0.0;
            m2_sum[i] += m2[i];
            m2[i] = 0.0;

        }

        // for (int i = 0; i < top_b; ++i)
        // {
        //     if((hist[i] != 0)) 
        //     std::cout << "g[" << i << "]=" << g[i] << std::endl;
        //     // std::cout << "hist = " << hist[i] << std::endl;
        // }

        g[0] = 0.0;

        it_count++;

        ss.str("");
        ss << "test_g/DoS-L=" << L;
        boost::filesystem::create_directories(ss.str().c_str());

        ss << "/" << it_count << ".dat";
        // std::cout << "write to " << ss.str() << std::endl;

        test_g_f.open(ss.str().c_str());

        for(int i = 0; i < top_b; i++)  {
            if((i!=2*(top_b)-1) && (hist[i] != 0))    {
                test_g_f << std::fixed << std::setprecision(6) << i << "\t" << energy(i) << "\t" << g[i] << "\n";
                std::cout << "G[" << i << "] = " << g[i]  << "; H[" << i << "]=" << hist[i] << std::endl;
            }
        }

        test_g_f.close();

        ss.str("");
        ss << "test_g/DoS-L=" << L << "/temp";
        boost::filesystem::create_directories(ss.str().c_str());

        ss.str("");
        ss << "test_g/DoS-" << "L=" << L << "/temp/" << it_count << ".plot";
        graph_g_f.open(ss.str().c_str());

        ss.str("");
        ss << "test_g/DoS-L=" << L << "/graphs";
        boost::filesystem::create_directories(ss.str().c_str());

        ss.str("");
        ss << "test_g/DoS-" << "L=" << L << "/" << it_count << ".dat";

        graph_g_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"test_g/DoS-L=" << L << "/graphs/" << it_count << ".jpg\"\n" << \
                 "set grid x y\n" << \
                 "set xlabel \"i\"\n" << \
                 "set ylabel \"G(i)\"\n" << \
                 "plot \"" << ss.str() << "\" using 1:3 title \"landau-wang-" << L << "-iteration-" << it_count << "\" with lines lt rgb \"red\"";

        graph_g_f.close();

        std::cout << std::fixed << std::setprecision(8) << "L - " << L << ": f = " << f << ", ln_f = " << ln_f << ", mcs = " << mcs << std::endl;

        f = pow(f, 0.5);    // Изменяем множитель

    }   

    ss.str("");
    ss << "plot_test_g_graph-L=" << L << ".sh";

    graph_sh_g.open(ss.str().c_str());
    graph_sh_g << "#!/bin/bash\n" <<  "gnuplot test_g/DoS-L=" << L << "/temp/*.plot\n";
    graph_sh_g <<  "convert -delay 100 -loop 0 test_g/DoS-L=" << L <<  \
                    "/graphs/{1..20}.jpg animate-DoS-L=" << L << ".gif\n";
    graph_sh_g.close();

    double EE, EE2, GE, Ut, Ft, St, Ct, Xt, MM, MM2, Mt, lambdatemp, lambda;
    double Ct_last, Xt_last, T_Ct_max, T_Xt_max;

    // out_f_mm_mm2 - файл для вывода намагниченности и восприимчивости
    std::ofstream out_f_td, out_f_ds, out_f_mm_mm2, out_t_max, plot_f, script_f, time_f;

    char * filename_out_td = new char [100];
    char * filename_out_ds = new char [100];
    char * filename_out_f_mm_mm2 = new char [100];
    char * filename_out_t_max = new char [100];

    // ---------------------------------------
    ss.str("");
    ss << "results/TermodinamicalStat_L=" << L << ".dat";

    strcpy(filename_out_td ,ss.str().c_str());

    out_f_td.open(filename_out_td);
    if(!out_f_td) std::cout << "Cannot open " << filename_out_td << ". Check permissions or free space!";
    out_f_td << "T\tUt\tFt\tSt\tCt\n";
    // ---------------------------------------
    ss.str("");
    ss << "results/DensityStat_L=" << L << ".dat";

    strcpy(filename_out_ds, ss.str().c_str());

    out_f_ds.open(filename_out_ds, std::ios::out);
    if(!out_f_ds) std::cout << "Cannot open " << filename_out_ds << ". Check permissions or free space!";
    out_f_ds << "i\tE(i)\tg[i]\thist[i]\n";
    // ---------------------------------------
    ss.str("");
    ss << "results/MagnetStats_L=" << L << ".dat";
    strcpy(filename_out_f_mm_mm2, ss.str().c_str());

    out_f_mm_mm2.open(filename_out_f_mm_mm2);
    if(!out_f_mm_mm2) std::cout << "Cannot open " << filename_out_f_mm_mm2 << ". Check permissions or free space!";
    // ---------------------------------------
    ss.str("");
    ss << "results/T_Critical_L=" << L << ".dat";
    
    strcpy(filename_out_t_max, ss.str().c_str());
    out_t_max.open(filename_out_t_max);
    if(!out_t_max) std::cout << "Cannot open " << filename_out_t_max << ". Check permissions or free space!";
    // ---------------------------------------
    for(int i = 0; i < top_b; i++)  {

        if(hist[i] != 0)  {

            m_sum[i]  /= it_count;
            m2_sum[i] /= it_count;

            // std::cout << "m_sum[" << i << "]=" << m_sum[i] << ";\tIT = " << it_count << ";" << std::endl;

        }

    }

    Ct_last = 0;
    Xt_last = 0;

    for(double T = 0.01; T <= 8; T += 0.01)  {

        EE = 0;
        EE2 = 0;
        GE = 0;
        MM = 0;
        MM2 = 0;

        lambda = 0;
        lambdatemp = 0;
        
        for(int i = 0; i < top_b; i++)  {
            if((i!=0) && i!=L*L-1 && hist[i]!=0)    {
                lambdatemp = g[i] - energy(i)/T;
                if(lambdatemp > lambda) lambda = lambdatemp;
            }
        }

        for(int i = 0; i < top_b; i++) {
            if((i!=1) && (i!=L*L-1) && (hist[i]!=0))    {

                EE  += energy(i)*exp(g[i]-(energy(i))/T-lambda);
                EE2 += energy(i)*energy(i)*exp(g[i]-(energy(i))/T-lambda);
                GE  += exp(g[i]-energy(i)/T-lambda);
                
                MM  += m_sum[i]*exp(g[i]-(energy(i))/T-lambda);

                // MM_mult_MM += m[i]*exp(g[i]-(energy(i))/T-lambda);

                // MM2 += m[i]*m[i]*exp(g[i]-(energy(i))/T-lambda);
                MM2 += m2_sum[i]*exp(g[i]-(energy(i))/T-lambda);

                // MM2 += m2[i];

            }
        }

        // MM2  = MM2/GE;

        // Ut = EE/GE;
        // Ft = -T*lambda-(T)*log(GE);
        // St = (Ut-Ft)/T;
        // Ct = ((EE2/GE)-Ut*Ut)/(T*T);

        // Mt = MM/GE;
        // Xt = ((MM2/GE)-Mt*Mt)/T;

        MM2 = MM2/GE;

        Ut = EE/GE;
        Ft = -T*lambda-(T)*log(GE);
        St = (Ut-Ft)/T;
        Ct = ((EE2/GE)-Ut*Ut)/(T*T);

        Mt = MM/GE;
        Xt = (Mt*Mt-(MM2/GE))/T;

        if(Xt > Xt_last)    {

            T_Xt_max = T;
            Xt_last = Xt;

        }

        if(Ct > Ct_last)    {

            T_Ct_max = T;
            Ct_last = Ct;

        }
        
        out_f_mm_mm2 << std::setprecision(3) << (float)T << std::setprecision(6) << "\t" << Mt/(L*L) << "\t" << Xt/(L*L) << std::endl; 
        // out_f_mm_mm2 << std::setprecision(3) << (float)T << std::setprecision(6) << "\t" << MM/(L*L) << "\t" << GE << std::endl; 

        // if((Ut==Ut)==true&&(Ft==Ft)==true&&(St==St)==true&&(Ct==Ct)==true) // Не пишем NaN 
        out_f_td << std::fixed << std::setprecision(4) << T << "\t" << Ut/(L*L) << "\t" << Ft/(L*L) << "\t" << St/(L*L) << "\t" << Ct/(L*L) << "\n";

    }

    for(int i = 0; i < top_b; i++)  {
        if((i!=2*(L*L)-1) && (hist[i] != 0))    {
            out_f_ds << std::fixed << std::setprecision(6) << i << "\t" << energy(i) << "\t" << g[i] << "\t" << hist[i] << "\n";
        }
    }

    std::cout << "T(Ct)_max = " << T_Ct_max << std::endl << \
                 "T(Xt)_max = " << T_Xt_max << std::endl;

    out_t_max << "T(Ct)_max = " << T_Ct_max << std::endl << \
                 "T(Xt)_max = " << T_Xt_max << std::endl;    

    ss.str("");
    ss << "temp/TermodinamicalStat_L=" << L;

    boost::filesystem::create_directories(ss.str().c_str());

    ss.str("");
    ss << "temp/TermodinamicalStat_L=" << L << "/Ut.plot";
    plot_f.open(ss.str().c_str());

    plot_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"graph/TermodinamicalStat_L=" << L << "/Ut.jpg\"\n" << \
                 "set grid x y\n" << \
                 "set xlabel \"T\"\n" << \
                 "set ylabel \"Ut\"\n" << \
                 "plot \"" << filename_out_td << "\" using 1:2 title \"landau-wang-" << L << "\" with lines lt rgb \"red\"";

    plot_f.close();

    ss.str("");
    ss << "temp/TermodinamicalStat_L=" << L << "/Ft.plot";
    plot_f.open(ss.str().c_str());

    plot_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"graph/TermodinamicalStat_L=" << L << "/Ft.jpg\"\n" << \
                 "set grid x y\n" << \
                 "set xlabel \"T\"\n" << \
                 "set ylabel \"Ft\"\n" << \
                 "plot \"" << filename_out_td << "\" using 1:3 title \"landau-wang-" << L << "\" with lines lt rgb \"red\"";

    plot_f.close();

    ss.str("");
    ss << "temp/TermodinamicalStat_L=" << L << "/St.plot";
    plot_f.open(ss.str().c_str());

    plot_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"graph/TermodinamicalStat_L=" << L << "/St.jpg\"\n" << \
                 "set grid x y\n" << \
                 "set xlabel \"T\"\n" << \
                 "set ylabel \"St\"\n" << \
                 "plot \"" << filename_out_td << "\" using 1:4 title \"landau-wang-" << L << "\" with lines lt rgb \"red\"";

    plot_f.close();

    ss.str("");
    ss << "temp/TermodinamicalStat_L=" << L << "/Ct.plot";
    plot_f.open(ss.str().c_str());

    plot_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"graph/TermodinamicalStat_L=" << L << "/Ct.jpg\"\n" << \
                 "set grid x y\n" << \
                 "set xlabel \"T\"\n" << \
                 "set ylabel \"Ct\"\n" << \
                 "plot \"" << filename_out_td << "\" using 1:5 title \"landau-wang-" << L << "\" with lines lt rgb \"red\"";

    plot_f.close();

    ss.str("");
    ss << "temp/DensityStat_L=" << L;

    boost::filesystem::create_directories(ss.str().c_str());

    ss.str("");
    ss << "temp/DensityStat_L=" << L << "/Ei.plot";
    plot_f.open(ss.str().c_str());

    plot_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"graph/DensityStat_L=" << L << "/Ei.jpg\"\n" << \
                 "set grid x y\n" << \
                 "set xlabel \"i\"\n" << \
                 "set ylabel \"E(i)\"\n" << \
                 "plot \"" << filename_out_ds << "\" using 1:2 title \"landau-wang-" << L << "\" with lines lt rgb \"red\"";

    plot_f.close();   

    ss.str("");
    ss << "temp/DensityStat_L=" << L << "/gi.plot";
    plot_f.open(ss.str().c_str());

    plot_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"graph/DensityStat_L=" << L << "/gi.jpg\"\n" << \
                 "set grid x y\n" << \
                 "set xlabel \"i\"\n" << \
                 "set ylabel \"g(i)\"\n" << \
                 "plot \"" << filename_out_ds << "\" using 1:3 title \"landau-wang-" << L << "\" with lines lt rgb \"red\"";

    plot_f.close();  

    ss.str("");
    ss << "temp/DensityStat_L=" << L << "/Hi.plot";
    plot_f.open(ss.str().c_str());

    plot_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"graph/DensityStat_L=" << L << "/Hi.jpg\"\n" << \
                 "set grid x y\n" << \
                 "set yrange [0:10000000]\n" << \
                 "set xlabel \"i\"\n" << \
                 "set ylabel \"H(i)\"\n" << \
                 "plot \"" << filename_out_ds << "\" using 1:4 title \"landau-wang-" << L << "\" with lines lt rgb \"red\"";

    plot_f.close();  

    ss.str("");
    ss << "graph/TermodinamicalStat_L=" << L;

    boost::filesystem::create_directories(ss.str().c_str());

    ss.str("");
    ss << "graph/DensityStat_L=" << L;

    boost::filesystem::create_directories(ss.str().c_str());

    ss.str("");
    ss << "plot_graph_L=" << L << ".sh";

    script_f.open(ss.str().c_str());

    ss.str("");
    ss << "#!/bin/bash\n";
    ss << "gnuplot temp/TermodinamicalStat_L=" << L << "/*.plot\n";
    ss << "gnuplot temp/DensityStat_L=" << L << "/*.plot\n";

    script_f << ss.str();

    // ss.str("");
    // ss << "sh plot_graph_L=" << L << ".sh";

    // chdir(".");
    // system(ss.str().c_str());

    time_e = omp_get_wtime();

    std::cout << "Time: " << time_e - time_b << "'s" << std::endl;

    ss.str("");
    ss << "Time-" << L << "-PP-1";

    time_f.open(ss.str().c_str());
    time_f << time_e - time_b << "'s" << std::endl;
    time_f.close();

    out_f_td.close();
    out_f_ds.close();
    out_f_mm_mm2.close();
    out_t_max.close();

    delete [] filename_out_td;
    delete [] filename_out_ds;
    delete [] filename_out_f_mm_mm2;
    delete [] filename_out_t_max;

    for(int i = 0; i < L; i++)  {
        delete spin[i];
    };
    delete [] spin;

    delete [] hist;
    delete [] g;
    delete [] m;
    delete [] m_sum;
    delete [] m2;
    delete [] m2_sum;

    return 0;

}

int LandauWangSerialHeisenberg()  {

    int i, j, k, N, a, b=0, b_new,
    n, skip, min_steps, mcs, rng, top_b, *hist;

    double f, fmin, dec_pow, flat_thres,
    r, ksi1, ksi2, ksi3, p,
    *g, *m, ***s_x, ***s_y, ***s_z,
    Delta[36]={0.01,
        0.522701295579425,
        0.581317791984826,
        0.63601096704108,
        0.686780820748186,
        0.733627353106146,
        0.776550564114959,
        0.815550453774625,
        0.850627022085143,
        0.881780269046515,
        0.90901019465874,
        0.932316798921817,
        0.951700081835748,
        0.967160043400532,
        0.978696683616168,
        0.986310002482658,
        0.99,
        0.989766676168195,
        0.985610030987244,
        0.977530064457145,
        0.9655267765779,
        0.949600167349507,
        0.929750236771967,
        0.90597698484528,
        0.878280411569447,
        0.846660516944466,
        0.811117300970338,
        0.771650763647063,
        0.728260904974641,
        0.680947724953072,
        0.629711223582357,
        0.574551400862494,
        0.515468256793484,
        0.452461791375326,
        0.385532004608023,
        0.314678896491571},
    arcsin_Theta=1.0;

    FILE *OUT;

    OUT = fopen("config.cfg", "r");

    if(fscanf(OUT, "%d", &Lx))  {
        printf("Cool! We can read Lx\n");
    };
    if(fscanf(OUT, "%d", &Ly))  {
        printf("Cool! We can read Ly\n");
    }
    if(fscanf(OUT, "%d", &Lz))  {
        printf("Cool! We can read Lz\n");
    }
    if(fscanf(OUT, "%d", &rng)) {
        printf("Cool! We can read rng\n");
    }
    if(fscanf(OUT, "%lf", &f))  {
        printf("Cool!We can read f\n");
    }
    if(fscanf(OUT, "%lf", &fmin))   {
        printf("Cool!We can read fmin\n");
    }
    if(fscanf(OUT, "%lf", &dec_pow))    {
        printf("Cool!We can read dec_pow\n");
    }
    if(fscanf(OUT, "%lf", &flat_thres)) {
        printf("Cool!We can read flat_thres\n");
    }
    if(fscanf(OUT, "%d", &min_steps))   {
        printf("Cool!We can read min_steps\n");
    }
    if(fscanf(OUT, "%d", &skip))    {
        printf("Cool!We can read skip\n");
    }

    fclose(OUT);

    // Lx = 8;
    // Ly = 8;
    // Lz = 8;

    // rng = 10000;

    // f = 2.718281828;
    // fmin = 1.00000001;

    // dec_pow = 0.5;
    // flat_thres = 0.8;

    // min_steps = 10000;
    // skip = 10000;

    N = Lx*Ly*Lz;

    /** UNCOMMENT FOR ANISOTROPY **/
    //arcsin_Theta = asin(1.0-Delta[Lz]);

    if(!(Lz%2))
        top_b=Lx*Ly*Lz+1;
    else
        top_b=Lx*Ly*(Lz-1)+1;//if that happens then everything is bad

    hist = new int [top_b];
    g = new double [top_b];
    m = new double [top_b];

    s_x = new double** [Lx];
    s_y = new double** [Lx];
    s_z = new double** [Lx];

    for(i=0;i<Lx;i++)
        {
            s_x[i] = new double* [Ly];
            s_y[i] = new double* [Ly];
            s_z[i] = new double* [Ly];

                for(j=0;j<Ly;j++)
                    {
                        s_x[i][j] = new double [Lz];
                        s_y[i][j] = new double [Lz];
                        s_z[i][j] = new double [Lz];
                    }
        }

    for(i=0;i<Lx;i++)
        {
            for(j=0;j<Ly;j++)
                {
                    for(k=0;k<Lz;k++)
                        {
                            s_x[i][j][k] = 0.0;
                            s_y[i][j][k] = 0.0;
                            s_z[i][j][k] = 1.0;
                        }
                }
        }

    for(i=0;i<top_b;i++)
        g[i]=1.0;

    CRandomMersenne Mersenne(time(0));

    while(f>fmin)
        {
            printf("f=%.9lf\n", f);

            double lnf=log(f);
            int c, cont=1;

            for(a=0;a<top_b;a++)
                hist[a] = 0;

            n=0;
            c=skip+1;
            mcs=0;

            do
                {
                    if(mcs%100000==0)
                        printf("mcs %d\n", mcs);

                    for(a=0;a<N;a++)
                        {
                            i = Mersenne.IRandomX(0, Lx-1);
                            j = Mersenne.IRandomX(0, Ly-1);
                            k = Mersenne.IRandomX(0, Lz-1);

                            b_new=b+(int)((s_x[i][j][k] * NeighboursHeisenberg(s_x,i,j,k)) +
                            (s_y[i][j][k] * NeighboursHeisenberg(s_y,i,j,k)) +
                            (s_z[i][j][k] * NeighboursHeisenberg(s_z,i,j,k)));

                           if(b_new >= 0 && b_new < top_b )
                                {
                                    p = exp(g[b] - g[b_new]);

                                    if((p>=1.0)||(Mersenne.Random() < p))
                                        {
                                            do
                                                {
                                                    ksi1 = (double) Mersenne.IRandomX(-rng,rng)/rng;
                                                    ksi2 = (double) Mersenne.IRandomX(-rng,rng)/rng;
                                                    ksi3 = (double) Mersenne.IRandomX(-rng,rng)/rng;

                                                    r = sqrt(ksi1*ksi1 + ksi2*ksi2 + ksi3*ksi3);
                                                } while(r>1.0);

                                            s_x[i][j][k] = ksi1*arcsin_Theta/r;
                                            s_y[i][j][k] = ksi2*arcsin_Theta/r;
                                            s_z[i][j][k] = ksi3/r;

                                            b = b_new;
                                        }
                                    g[b]+=lnf;
                                    m[b]+=(s_x[i][j][k]*s_y[i][j][k]*s_z[i][j][k]);

                                    hist[b] += 1;
                                    n++;
                                }
                        }
                    mcs++;
                    c++;

                    if((mcs>=min_steps) && (c>=skip))
                        {
                            int div;
                            c=0;
                            div=top_b>Lx*Ly*Lz-1 ? top_b-2 : top_b-1;

                            for(a=0;
                                (a<top_b) && ((a==1) || (a==Lx*Lx*Lz-1) || ((double)hist[a]/n*div > flat_thres));
                                    a++);
                            if(a==top_b)
                                cont=0;
                            else
                                cont=1;
                        }
                    else
                        cont=1;

                } while(cont);
                /** alternative flatness check **/
                    /*if((mcs >= min_steps) && (c >= skip))
                        {
                            c = 0;

                            int h_count = 0;
                            double h_delt;
                            double h_sum = 0.0;

                            for(int i = 0; i < top_b; i++)
                                {
                                    if(hist[i] != 0)
                                        {
                                            h_sum += (double)hist[i]/min_steps;
                                            h_count++;
                                        }

                                }

                            cont = 0;

                            for (int i = 0; i < top_b; i++)
                                {
                                    if(hist[i] != 0)
                                        {
                                            h_delt = (double)hist[i]/(h_sum/h_count)/min_steps;
                                            if((h_delt <= flat_thres) || (h_delt >= 1.0+flat_thres))
                                                cont = 1;
                                        }
                                }
                        }
                    else
                        cont = 1;
                }while(cont);*/

        for(a=1;a<top_b;a++){
            g[a] -= g[0];
            m[i] = abs(m[i]);
        }

        g[0] = 0.0;

        f = pow(f,dec_pow);
        }

    OUT = fopen("data.dat", "w+");
        {
            for(a=0;a<top_b;a++)
                {
                    if((a!=1)&&(a!=Lx*Ly*Lz-1))
                        {
                            fprintf(OUT, "%d", a-top_b/2);
                            fprintf(OUT, "\t%.9lf", g[a]+log(100000.0));
                            fprintf(OUT, "\t%d\n", hist[a]);
                            fprintf(OUT, "\t%.9lf\n", m[a]);
                        }
                }
        }
    fclose(OUT);

    OUT = fopen("temp.cfg","w+");
    fprintf(OUT, "%d\n", Lx);
    fprintf(OUT, "%d\n", Ly);
    fprintf(OUT, "%d\n", Lz);
    fprintf(OUT, "%d\n", top_b);
    fclose(OUT);

    printf("\a \a");

    return 0;

}

int DataProcessing()    {


    int T_step_max = 800;
    double T0 = 0.1, T1 = 8.1;

    int N, n;
    int i, T_step;
    double prom, *ln_g, *H, *E;
    double E_eq, T, dT, Z, E2_eq;
    char fname[80];
    FILE *OUT;

    OUT = fopen("temp.cfg","r");
    if(fscanf(OUT, "%d", &Lx))  {
        printf("Cool! We can read Lx\n");
    };
    if(fscanf(OUT, "%d", &Ly))  {
        printf("Cool! We can read Ly\n");
    };
    if(fscanf(OUT, "%d", &Lz))  {
        printf("Cool! we can read Lz\n");
    };
    if(fscanf(OUT, "%d", &n))   {
        printf("Cool! we can read \"n = %d\"\n", n);
    };
    fclose(OUT);

    N = Lx*Ly*Lz;

    ln_g = new double [n];
    H = new double [n];
    E = new double [n];

    dT = (T1-T0)/T_step_max;

    OUT = fopen("data.dat","r");

    for(i=0;i<n;++i)
    {
        if(fscanf(OUT,"%lf", &E[i]))    {
            printf("E[%d] = %.8lf\t", i, E[i]);
        };
        if(fscanf(OUT,"%lf", &ln_g[i])) {
            printf("ln_g[%d] = %.8lf\t", i, ln_g[i]);
        };
        if(fscanf(OUT,"%lf", &H[i]))    {
            printf("H[%d] = %.8lf\n", i, H[i]);
        };
    }

    fclose(OUT);

    OUT = fopen("result.dat","w+");
    for(T_step=0;T_step<=T_step_max;++T_step)
    {
        T = T0 + dT*T_step;
        printf("%.4lf", T);
        fprintf(OUT,"%.4lf", T);

        prom = ln_g[0] - E[0]/T;

        for(i=1;i<n;++i)
        {
            if((ln_g[i] - E[i]/T)>prom)
            {
                prom = ln_g[i] - E[i]/T;
            }
        }

        if(T_step%10==0)
        {
            printf("T = %.5lf", T);
            printf("\n");
        }

        E_eq = 0.0;
        E2_eq = 0.0;
        Z = 0.0;

        for(i=0;i<n;++i)
        {
            E_eq = E_eq + E[i]*exp(ln_g[i]-E[i]/T - prom);
            E2_eq = E2_eq + E[i]*E[i]*exp(ln_g[i]-E[i]/T - prom);
            Z = Z + exp(ln_g[i]-E[i]/T - prom);
        }
        fprintf(OUT, "\t%.9lf", (E_eq/Z)/N);
        fprintf(OUT, "\t%.9lf", (((E2_eq/Z)-(E_eq/Z)*(E_eq/Z))/(T*T))/N);
        fprintf(OUT, "\t%.9lf", (-T*prom-(T)*log(Z))/N);
        fprintf(OUT, "\t%.9lf", (((E_eq/Z)-(-T*prom-(T)*log(Z)))/T)/N);
        fprintf(OUT, "\n");
    }
    fclose(OUT);

    delete(H);
    delete(ln_g);
    delete(E);
    return 0;

}

double NeighboursHeisenberg( double ***s, int i, int j, int k )
{
    double res;

    if(i==0) res = s[Lx-1][j][k];
    else res = s[i-1][j][k];

    if(i==Lx-1) res = res + s[0][j][k];
    else res = res + s[i+1][j][k];

    if(j==0) res = res + s[i][Ly-1][k];
    else res = res + s[i][j-1][k];

    if(j==Ly-1) res = res + s[i][0][k];
    else res = res + s[i][j+1][k];

    if(k==0) res = res + s[i][j][Lz-1];
    else res = res + s[i][j][k-1];

    if(k==Lz-1) res = res + s[i][j][0];
    else res = res + s[i][j][k+1];

    return res;
}

inline int neighbour_spins(int i, int j)    {

    int result;

    if(i==0)    result=spin[L-1][j];    
    else        result=spin[i-1][j];
    if(i==L-1)  result+=spin[0][j];
    else        result+=spin[i+1][j];
    if(j==0)    result+=spin[i][L-1];
    else        result+=spin[i][j-1];
    if(j==L-1)  result+=spin[i][0];
    else        result+=spin[i][j+1];

    return result;
}
