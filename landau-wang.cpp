/*
 * landau-wang.cpp
 *
 *  Created on: 24 мая 2014 г.
 *      Author: nlare
 */

#include <iostream>
#include <iomanip>
#include <cstdlib>
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

int neighbour_spins(int,int);

int main(int argc, char *argv[])  {

    int *hist;              // Массив гистограммы энергий, т.е количества посещений данного энергетического состояния
    double *g;              /* Массив плотности состояний
                             * Изначально каждый элемент массива принимается равным единице
                             */

	double *m;              /* Массив для намагниченности
                             * Изначально каждый элемент массива принимается равным нулю
                             */
	
	double *m2;				/* Массив для квадрата намагниченности
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

    L = 16;
    f = 2.7182818284;  // В работе Ландау-Ванга было указано значение "f" равное экспоненте 
    f_min = 1.000001;  // Данная переменная должна быть около единицы
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
    m2 = new double [4*L*L];
    
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
    }

	for(int i = 0; i < 2*L*L; i++)    {
        m2[i] = 0.0;
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
                    
                    g[b] += ln_f;   	  // Увеличиваем текущий энергетический уровень на логарифм "f"
					m[b] += spin[ci][cj]; // Увеличиваем значение намагниченности для данного энергетического уровня

					// m2[b] += spin[ci][cj]*spin[ci][cj];   // Квадрат намагниченности для расчета восприимчивости

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

            m[i] = abs(m[i]);

            m2[i] = abs(m2[i]);

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

    // out_f_mm_mm2 - файл для вывода намагниченности и восприимчивости
    std::ofstream out_f_td, out_f_ds, out_f_mm_mm2, plot_f, script_f, time_f;

    char * filename_out_td = new char [100];
    char * filename_out_ds = new char [100];
    char * filename_out_f_mm_mm2 = new char [100];

    ss.str("");
    ss << "results/TermodinamicalStat_L=" << L << ".dat";

    strcpy(filename_out_td ,ss.str().c_str());

    out_f_td.open(filename_out_td);
    if(!out_f_td) std::cout << "Cannot open " << filename_out_td << ".Check permissions or free space";
    out_f_td << "T\tUt\tFt\tSt\tCt\n";

    ss.str("");
    ss << "results/DensityStat_L=" << L << ".dat";

    strcpy(filename_out_ds, ss.str().c_str());

    out_f_ds.open(filename_out_ds, std::ios::out);
    if(!out_f_ds) std::cout << "Cannot open " << filename_out_ds << ".Check permissions or free space";
    out_f_ds << "i\tE(i)\tg[i]\thist[i]\n";

    ss.str("");
	ss << "results/MagnetStats_L=" << L << ".dat";
	strcpy(filename_out_f_mm_mm2, ss.str().c_str());

	out_f_mm_mm2.open(filename_out_f_mm_mm2);

    double m_min = abs(m[0]);
    double m2_min = abs(m2[0]);

  // //       std::cout << "m[0] = " << m_min << std::endl;

        for(int i = 0; i < top_b; i++)  {
          
          if(hist[i] != 0)  {

            m[i] = abs(m[i]);

            m2[i] = abs(m2[i]);

            if((m[i] > 0) && (m[i] < m_min)) m_min = m[i];

            if((m[i] > 0) && (m2[i] < m2_min)) m2_min = m2[i];

             // std::cout << "m_min = " << m_min << "; m[" << i << "]=" << m[i] << "; m2[" << i << "]=" << m2[i] << std::endl;
          }

        }

        for(int i = 0; i < top_b; i++)  {

            if(hist[i] != 0)  {

              m[i]  -= m_min;
              m2[i] -= m2_min;

  //             // std::cout << "m_min = " << m_min << "; m[" << i << "]=" << m[i] << "; m2[" << i << "]=" << m2[i] << std::endl;

            }

        }

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
				
				MM  += m[i]*exp(g[i]-(energy(i))/T-lambda);

                // MM_mult_MM += m[i]*exp(g[i]-(energy(i))/T-lambda);

				MM2 += m[i]*m[i]*exp(g[i]-(energy(i))/T-lambda);

                // MM2 += m2[i];

            }
        }

        // MM_mult_MM = MM_mult_MM/GE;

		// MM2	= MM2/GE;

        Ut = EE/GE;
        Ft = -T*lambda-(T)*log(GE);
        St = (Ut-Ft)/T;
        Ct = ((EE2/GE)-Ut*Ut)/(T*T);

        Mt = MM/GE;
        Xt = ((MM2/GE)-Mt*Mt)/T;
		
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

    delete [] filename_out_td;
    delete [] filename_out_ds;
	delete [] filename_out_f_mm_mm2;

    for(int i = 0; i < L; i++)  {
        delete spin[i];
    };
    delete [] spin;
    delete [] hist;
    delete [] g;
    delete [] m;
    delete [] m2;
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
