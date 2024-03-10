using namespace std;

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <random>
#include <chrono>


/* Ponizsze parametry zmieniac zaleznie od potrzeb */
//liczebnosc populacji
#define POPSIZE 25
//maksymalna liczba pokolen
#define MAXGENS 250
//liczba zmiennych w zadaniu
#define NVARS 30
//prawdopodobienstwo krzyzowania
#define PXOVER 0.8
//prawdopodobienstwo mutacji
#define PMUTATION 0.1
//rzeczywista liczba genow
#define NBIN 48
//delta rozkladu normalnego
#define delta 10
// liczba PI
#define M_PI 3.14159

//wartosc maksymalna (stala C)
#define MAX 25

//mean
auto mean=0;
//stddev
auto stddev = 0.2;

#define TRUE 1
#define FALSE 0

//genotyp, czlonek populacji
 struct genotype
 {
  double gene[NVARS];			    //lancuch zmiennych
  double fitness;					//dopasowanie genotypu
  double upper[NVARS]; 		        //gorne ograniczenie na zmienne genotypu
  double lower[NVARS];		        //dolne ograniczenie na zmienne genotypu
  double rfitness;				    //dopasowanie wzgledne
  double cfitness;				    //dopasowanie laczne
  double genebin[NBIN];             //wartsci binarne genow
  double oblicz_dopasowanie();      //metoda obliczajaca dopasowanie
 };

//numer biezacego pokolenia
int generation;
//najlepszy osobnik
int cur_best;
//plik wyjsciowy
FILE *galog;

//populacja
genotype population[POPSIZE+1];
//nowa populacja zamieniajaca stare pokolenie
genotype newpopulation[POPSIZE+1];
//srednia funkcja przystosowania poprzedniego pokolenia
double	avg_pop;
//licznik zegara
double zegar;
//ograniczenia dolne i gorne zakresu zmiennych
double lbound, ubound;

/**************************************************************************/
// funkcja optymalizowana
/**************************************************************************/

double genotype::oblicz_dopasowanie()
{
    int k;
    double w1, w2;
    double sum, sumcos;
    double c = 2 * 3.14159265358979323846;

	double dopasowanie = 0;

    sum=0;
    sumcos=0;
    for (k=0; k<NVARS; k++)
    {
        sum = sum + gene[k] * gene[k];
        sumcos = sumcos + cos(c * gene[k]);
    }
    w1 = -20 * exp(-0.2 * sqrt(sum/NVARS));
    w2 = exp(sumcos/NVARS);

    //dopasowanie = 100*(gene[1]-gene[0]*gene[0])*(gene[1]-gene[0]*gene[0]) + (gene[0]-1)*(gene[0]-1);
    dopasowanie = w1 - w2 + 20 + exp(1);
return dopasowanie;
};

/**************************************************************************/
//Wymaga C++11
//Settings/Compiler/Have g++ follow C++11
/**************************************************************************/
double losuj(double m, double s)
{
    double wynik;
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();

    default_random_engine generator(seed);
    normal_distribution<double> distribution(m, s);

    wynik = distribution(generator);
    //cout << "losowa " << wynik << endl;
    return wynik;
}

/*************************************************************************/
/*Procedura inicjujaca, nadaje genom wartosci wewnatrz ograniczen.       */
/*Ustawia na zero wartosci dopasowan dla wszystkich czlonkow popu-      */
/*lacji. Wczytuje ograniczenia dolne i gorne zmiennych z pliku           */
//wejsciowego "gadata.txt", nastepnie generuje losowo wartosci           */
//wewnatrz tych ograniczen dla kazdego genu w kazdym genotypie w         */
//populacji. Format pliku wejsciowego "gadata.txt" jak ponizej:          */
/*ograniczenie_dolne_zmiennej1  ograniczenie_gorne_zmiennej1             */
/*ograniczenie_dolne_zmiennej2  ograniczenie_gorne_zmiennej2 itd...      */
/*************************************************************************/

void initialize(void)
{
 FILE* infile;
 FILE* popfile;
 int i, j;
 double tmp;
// double temp[NVARS+1];
// int li;
// int bin[48];
// double liczba;


 if ((infile=fopen("gadata.txt","r"))==NULL)
    {
     fprintf(galog,"\nNie moge otworzyc pliku wejsciowego!/n");
     exit(1);
    }
 if ((popfile=fopen("popdata.txt","r"))==NULL)
    {
     fprintf(galog,"\nNie moge otworzyc pliku!/n");
     exit(2);
    }
 //ustalanie wartosci wewnatrz ograniczen
//fscanf(infile, "%lf", &lbound);
//fscanf(infile, "%lf", &ubound);
 for (i=0; i<NVARS; i++)
    {
     fscanf(infile, "%lf", &lbound);
     fscanf(infile, "%lf", &ubound);

     for (j=0; j<POPSIZE; j++)
        {
         population[j].fitness=0;
         population[j].rfitness=0;
         population[j].cfitness=0;
         population[j].lower[i]=lbound;
         population[j].upper[i]=ubound;
        }
    }
/*
for (i=0; i<NVARS; i++)
    {
     cout << "NVARS " << i << endl;
     for (j=0; j<POPSIZE; j++)
        {
         cout << population[j].lower[i] << " ";
         cout << population[j].upper[i] << endl;
//         population[j].gene[i]=rand()%2;
        }
    }
*/
/*
for (i=0; i<NVARS; i++)
  {
  fscanf(popfile, "%lf", &tmp);
  temp[i]=tmp;
  }
*/
/*
for (i=0; i<NVARS; i++)
  {
  cout << "temp " << temp[i] << " ";
  }
*/
for (i=0; i<POPSIZE; i++)
{
  for (j=0; j<NVARS; j++)
  {
  fscanf(popfile, "%lf", &tmp);
  population[i].gene[j]=tmp;
  }
}
 fclose(infile);
 fclose(popfile);

/*
cout << "po inicjalizacji " << endl;

for (j=0; j<POPSIZE; j++)
  {
    for (i=0; i<NVARS; i++)
    {
      cout << population[j].gene[i] << " ";
    }
    cout << endl;
  }
*/
//cout << "initialize OK" << endl;

}

/*************************************************************************/
/*Generator liczb losowych generujacy wartosci wewnatrz                  */
/*ograniczen                                                             */
/*************************************************************************/

 double randval(double low, double high)
 {
  double val;
  val=((double)(rand()%1000/1000.0)*(high - low)+low);
  return (val);
 }

/*************************************************************************/
/*Funkcja oceny ustalana przez uzytkownika.                              */
/*Po kazdej zmianie funkcji nalezy ponownie skompilowac program.         */
/*Obecnie jest to funkcja: x[1]^2-x[1]*x[2]+x[3]                         */
/*************************************************************************/

void evaluate (void)
{
 int mem;
 //int i;
 //double x[NVARS+1];

 for (mem=0; mem<POPSIZE; mem++)
    {
/*
     for (i=0; i<NVARS; i++)
        {
         x[i+1]=population[mem].gene[i];
        }

     for (i=0; i<NVARS; i++)
        {
         cout << x[i+1] << " ";
        }
     cout << "koniec iksow " << endl;
*/

     population[mem].fitness= MAX - population[mem].oblicz_dopasowanie();

     if (population[mem].fitness < 0) population[mem].fitness = 0;
    }
/*
cout << " fitness   " << endl;
for (i=0; i<POPSIZE; i++)
{
cout << population[i].fitness << "    " << population[i].gene[0] << " ";
cout << population[i].gene[1] << endl;
}
*/
//cout << endl << "evaluate OK" << endl;
}
/*************************************************************************/
/*Procedura keep_the_best. Zapamietuje ona najlepszego dotychcza-        */
/*wego czlonka populacji. Uwaga: w ostatnim elemencie macierzy           */
/*population znajduje sie kopia najlepszego osobnika.                    */
/*************************************************************************/

void keep_the_best (void)
{
 int mem;
 int i;
 cur_best =0;	//zapamietanie wskaznika najlepszego osobnika
 for (mem=0; mem<POPSIZE; mem++)
    {
     if (population[mem].fitness>population[POPSIZE].fitness)
       {
        cur_best=mem;
        population[POPSIZE].fitness=population[mem].fitness;
       }
    }
/*Po znalezieniu najlepszego czlonka populacji skopiuj geny              */

 for( i=0; i<NVARS; i++)
    {
     population[POPSIZE].gene[i]=population[cur_best].gene[i];
    }
/*
cout << "Skopiowane geny" << endl;
for (i=0; i<NVARS; i++)
  {
    cout << population[POPSIZE].gene[i] << " ";
  }
cout << endl;
*/
//cout << "keep th best OK" << endl;

}

/*************************************************************************/
/*Procedura elitist. Najlepszy osobnik z poprzedniego pokolenia          */
/*jest zapamietywany jako ostatni w macierzy. Jezeli najlepszy           */
/*osobnik z biezacego pokolenia jest gorszy niz najlepszy osobnik        */
/*z poprzednich pokolen, to ten ostatni zastepuje najgorszego            */
/*osobnika biezacej populacji                                            */
/*************************************************************************/

void elitist (void)
{
 int i;
 double best, worst;						//najlepsza i najgorsza wartosc dopasowania
 int best_mem, worst_mem;			    	//wskazniki do najlepszego i
 											//najgorszego osobnika
 for (i=0; i<POPSIZE-1; ++i)
    {
     if (population[i].fitness>population[i+1].fitness)
       {
        if (population[i].fitness>=best)
          {
           best=population[i].fitness;
           best_mem=i;
          }
        if (population[i+1].fitness<=worst);
          {
           worst=population[i+1].fitness;
           worst_mem=i+1;
          }
       }
      else
       {
        if (population[i].fitness<=worst);
          {
           worst=population[i].fitness;
           worst_mem=i;
          }
        if (population[i+1].fitness>=best);
          {
           best=population[i+1].fitness;
           best_mem=i+1;
          }
       }
    }

/*************************************************************************/
/*Jezeli najlepszy osobnik z nowej populacji jest lepszy niz             */
/*najlepszy osobnik z poprzednich populacji, to skopiuj                  */
/*najlepszego z nowej populacji, jezeli nie, to zastap                   */
/*najgorszego osobnika z biezacej populacji przez najlepszego            */
/*z poprzednich pokolen                                                  */
/*************************************************************************/

     if (best>=population[POPSIZE].fitness)
       {
        for (i=0; i<NVARS; i++)
           {
            population[POPSIZE].gene[i]=population[best_mem].gene[i];
           }
        population[POPSIZE].fitness=population[best_mem].fitness;
       }
     else
       {
        for (i=0; i<NVARS; i++)
           {
            population[worst_mem].gene[i]=population[POPSIZE].gene[i];
           }
        population[worst_mem].fitness=population[POPSIZE].fitness;
       }
/*
cout << endl;
for (i=0; i<NVARS; i++)
  {
    cout << population[POPSIZE].gene[i] << " ";
  }
*/
//cout << endl << "Elitist OK" << endl;

}

/*************************************************************************/
/*Procedura wyboru. Selekcja standardowa proporcjonalna do zadania       */
/*maksymalizacji realizujaca model elitarny - zapewnia, ze zawsze        */
/*przezywa najlepszy osobnik populacji.                                  */
/*************************************************************************/

void select (void)
{
 int mem, i, j;
 double sum=0;
 double p;
/*Wyznaczenie calkowitego dopasowania populacji                          */

 for (mem=0; mem<POPSIZE; mem++)
    {
     sum+=population[mem].fitness;
    }
/*Obliczenie dopasowania wzglednego.                                     */

 for (mem=0; mem<POPSIZE; mem++)
    {
     population[mem].rfitness=population[mem].fitness/sum;
    }
/*Obliczanie dopasowania lacznego.                                       */
population[0].cfitness=population[0].rfitness;
 for (mem=1; mem<POPSIZE; mem++)
    {
     population[mem].cfitness=population[mem-1].cfitness+
                                population[mem].rfitness;
    }
/*Koncowy wybor osobnikow przezywajacych na podstawie dopasowan lacznych */

 for( i=0; i<POPSIZE; i++)
    {
     p=rand()%1000/1000.0;
     if (p<population[0].cfitness)
       newpopulation[i]=population[0];
     else
       {
        for (j=0; j<POPSIZE; j++)
           {
            if (p>=population[j].cfitness&&p<population[j+1].cfitness)
              newpopulation[i]=population[j+1];
           }
       }
    }
/*Kopiowanie populacji po jej utworzeniu.                                */

 for (i=0; i<POPSIZE; i++)
    {
     population[i]=newpopulation[i];
    }
//cout << "Selection OK" << endl;
}

/*************************************************************************/
/*Funkcja pomocnicza wymieniajaca wzajem dwie zmienne.                   */
/*************************************************************************/

void swap (double *x, double *y)
{
 double temp;
 temp=*x;
 *x=*y;
 *y=temp;
}

/*************************************************************************/
/*Krzyzowanie. Wykonanie krzyzowania dwojga wybranych rodzicow.          */
/*************************************************************************/

void Xover (int one, int two)
{
 int i;
 int point;	//punkt krzyzowania
/*Wybor punktu krzyzowania                                               */

 if (NVARS>1)
   {
    if (NVARS == 2)
      {
       point = 1;
      }
    else
      {
       point = (rand() % (NVARS - 1) +1);
      }
/*Wykonanie krzyzowania.                                                 */

 for (i=0; i<point; i++)
    {
     swap(&population[one].gene[i], &population[two].gene[i]);
    }
  }
}


/*************************************************************************/
/*Wybor do krzyzowania. Wybor dwojga rodzicow, ktorzy wezma udzial       */
/*w krzyzowaniu. Funkcja realizuje krzyzowanie jednopunktowe.            */
/*************************************************************************/

void crossover (void)
{
 int mem, one;
 int first=0;	//obliczenie liczby wybranych osobnikow.
 double x;
 for (mem=0; mem<POPSIZE; ++mem)
    {
     x=rand()%1000/1000.0;
     if (x<PXOVER)
       {
        ++first;
        if (first % 2 == 0)
          {
           Xover (one, mem);
          }
        else
          {
           one=mem;
          }
       }
    }
 }

/*************************************************************************/
/*Mutacja. Losowa mutacja jednorodna. Zmienna wybrana do mutacji jest    */
/*zamieniana na wartosc losowa zawarta wewnatrz ograniczenia dolnego i   */
/*gornego tej zmiennej.                                                  */
/*************************************************************************/

void mutate (void)
{

 int i, j;
 double x;
// double liczba;
// int lp[NVARS][16];
// int bin[48];

//cout << "granice " << lbound << " " << ubound << endl;

 for (i=0; i<POPSIZE; i++)
    {
     for (j=0; j<NVARS; j++)
        {
         x=rand() % 1000/1000.0;
         if (x<PMUTATION)
           {
           //population[i].gene[j]=randval(lbound, ubound);

           population[i].gene[j]=population[i].gene[j] + losuj(mean, stddev);

           //cout << "nowy gen " << population[i].gene[j] << endl;
           /*
           population[i].gene[j]=population[i].gene[j]+gauss(delta);
           if (population[i].gene[j]<population[i].lower[j])
              {
              population[i].gene[j]=population[i].lower[j];
              }
           if (population[i].gene[j]>population[i].upper[j])
              {
              population[i].gene[j]=population[i].upper[j];
              }
           */
           }
        }
    }
//    cout << endl;
}
/*************************************************************************/
/*Funkcja obliczajaca i zapisujaca wyniki do pliku. Zapisywane sa        */
/*postepy symulacji. Dane zapisywane do pliku wyjsciowego sa oddzielane  */
/*przecinkami.                                                           */
/*************************************************************************/
/*
void report (void)
{

 int i;
 double best_val; 	//najlepsze dopasowanie populacji
 double avg;			//srednie dopasowanie populacji
 double stddev;		//odchylenie standardowe dopasowania populacji
 double sum_square;	//suma kwadratow do obliczania odchylenia standardowego
 double square_sum;	//kwadrat sumy do obliczania odchylenia standardowego
 double sum;			//calkowite dopasowanie populacji

//Obliczanie parametrow uzyskanych w wyniku symulacji
 sum=0.0;
 sum_square=0.0;
 for (i=0; i<POPSIZE; i++)
    {
     sum+=population[i].fitness;
     sum_square+=population[i].fitness*population[i].fitness;
    }
 avg=sum/(double)POPSIZE;
 square_sum=avg*avg*(double)POPSIZE;
 stddev=sqrt((sum_square-square_sum)/(POPSIZE-1));
 best_val=population[POPSIZE].fitness;
// fprintf(galog, "\n%5d, %6.3f, %6.3f, %6.3f, %6.3f \n\n",generation,best_val, avg, stddev/avg*100, (best_val-avg)/avg*100);
avg_pop=avg;
}
*/



/*************************************************************************/
/*Program glowny. W kazdym pokoleniu nastepuje wybor najlepszych         */
/*osobnikow, krzyzowanie i mutacja, a nastepnie ocena kolejnej populacji */
/*az do spelnienia warunku koncowego.                                    */
/*************************************************************************/

int main(int argc, char *argv[])
{
 int i;

  srand(time(0));

 if ((galog=fopen("galog.txt","a"))==NULL)
   {
    exit(3);
   }

 generation = 0;

 initialize();

 evaluate();

 keep_the_best();

 //while (population[POPSIZE].fitness<= MAX - 0.01)
 while (generation<1000000)
      {
       generation++;
       select();
       crossover();
       mutate();
//       report();
       evaluate();
//       keep_the_best();
       elitist();


   		//if (generation%100000==0) cout << generation  << " " << population[POPSIZE].fitness << endl;
        if (generation%10000==0) printf("%6.6f \n", MAX - population[POPSIZE].fitness);
    }

fprintf (galog, "\n\nkoniec = %6d ", generation);
cout << "koniec = " << generation << endl;


fprintf (galog, "\nSymulacja zakonczona \n");
fprintf (galog, "\nNajlepszy osobnik: \n");

for (i=0; i<NVARS; i++)
    {
    fprintf (galog,"\n var(%d)=%3.3f",i,population[POPSIZE].gene[i]);
    }


fprintf (galog,"\n\nNajlepsze dopasowanie=%3.3f",population[POPSIZE].fitness);

//Generowanie populacji poczatkowej
/*
for (i=0; i<POPSIZE; i++)
{
    fprintf (galog,"\n %3.6f %3.6f",randval(0, 3.14159), randval(0 , 3.14159));
}
*/
fclose (galog);

 for (i=0; i<NVARS; i++)
    {
     cout << population[POPSIZE].gene[i] << " ";
    }
 cout << endl;
 cout << population[POPSIZE].fitness << endl;
 zegar=clock()/1000000.0;
 cout << "Czas wykonania [s] " << zegar << endl;


 printf ("Sukces\n");

return EXIT_SUCCESS;
}




