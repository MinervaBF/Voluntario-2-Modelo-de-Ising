#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <fstream>
//Para usar la misma dimensión en todo el programa definimos dicha dimensión N.
#define N 32

using namespace std;

int main (void)
{

long int i, j, k, m1, m2, contador, icorrelacion, cont;
float temp, energ, p, num2, magnetizacion, energiamedia, energiacuadrado, caloresp, fcorrelacion, h;
int s[N+2][N+2];
double x[N+1][N+1];
//Los índices de la matriz van desde 0 a N+1 para las condiciones de contorno periódicas,
//si queremos una red de dimensión N, tendremos una matriz de N+2 elementos en cada fila.

FILE *f1, *f2, *f3, *f4;

//Abrimos los archivos
f1=fopen("magnetizacionmedia.txt","w");
f2=fopen("energiamedia.txt","w");
f3=fopen("calorespecifico.txt","w");
f4=fopen("tempfcorrelacion.txt","w");

srand(time(NULL));

//Introducimos el paso
h=(3.5-1.5)/9;

for(cont=0; cont<=9; cont++)
{

// Damos un valor a la temperatura
temp=1.5+cont*h;

// Generamos un estado inicial ordenado.

    for(i=0;i<=N;i++)
        for(j=0;j<=N;j++)
            s[i][j]=1;


//Forzamos que se cumplan las condiciones de contorno periódicas.

for(i=0;i<=N;i++)
{
    //Columnas
    s[0][i]=s[N][i];
    s[N+1][i]=s[1][i];
    //Filas
    s[i][0]=s[i][N];
    s[i][N+1]=s[i][1];
}


// Inicializo el vector en el que calculare el promedio de s(n,m)s(n+i,m) para la funcion de correlacion

for(i=0;i<=N;i++)
        for(j=0;j<=N;j++)
            x[i][j]=0.;


//Inicializo las variables.

magnetizacion=0.;
energiamedia=0.;
energiacuadrado=0.;
caloresp=0.;
fcorrelacion=0.;
contador=0;

// Para la función de correlación, damos un valor a la i que aparece en la fórmula.
icorrelacion=4;

// Comienzo el bucle que hará evolucionar el sistema.

for(k=1; k<=(10000*N*N); k=k+1) // Dejamos evolucionar el sistema 10^6 pMC
{


//Elegimos un punto al azar de la matriz. Para ello generamos dos números al azar dentro
//de las dimensiones de la matriz.
m1=1+(rand()%(N));
m2=1+(rand()%(N));

//Evaluamos p.
energ=2.*s[m1][m2]*(s[m1+1][m2]+s[m1-1][m2]+s[m1][m2+1]+s[m1][m2-1]);
p=exp(-1.*energ/temp);
if(1<p) p=1.;

//Generamos un número aleatorio perteneciente a [0,1] y si es menor que p cambiamos el signo del espín.

num2=(double)rand()/(double)RAND_MAX;
if (num2<p) s[m1][m2]=-1*s[m1][m2];

//Tenemos que volver a imponer que se cumplan las condiciones de contorno periódicas.
for(i=0;i<=N;i++)
{
    //Columnas
    s[0][i]=s[N][i];
    s[N+1][i]=s[1][i];
    //Filas
    s[i][0]=s[i][N];
    s[i][N+1]=s[i][1];
}

    // Para los cálculos, tomamos las medidas cada 100pMC.

   //if(k%(100*(N*N))==0)
    //{
        // Para la magnetización promedio:
        for (i=1; i<=N; i++)
        {
            for (j=1;j<=N; j++)
            {
                magnetizacion=magnetizacion+s[i][j];
                energiamedia=energiamedia-0.5*s[i][j]*(s[i][j+1]+s[i][j-1]+s[i+1][j]+s[i-1][j]);
                energiacuadrado=energiacuadrado+pow((0.5*s[i][j]*(s[i][j+1]+s[i][j-1]+s[i+1][j]+s[i-1][j])),2);
                if((i+icorrelacion) <= N)
                    x[i][j]=x[i][j]+s[i][j]*s[i+icorrelacion][j];
                else
                    x[i][j]=x[i][j]+s[i][j]*s[i+icorrelacion-N][j];
            }
        }
    contador=contador+1;
    //}


}

// Finalicemos los cálculos de las variables promediadas.

// Para la magnetización promedio:
magnetizacion=(1./(pow(N,2)*contador))*abs(magnetizacion);


// Para la energía media.
energiamedia=energiamedia/(2.*N*contador);

// Para el calor específico.
caloresp=(1./pow(N,2)*temp)*(energiacuadrado/contador-pow(energiamedia,2));

// Para la función de correlación.
for (i=1; i<=N; i++)
            for (j=1;j<=N; j++)
                x[i][j]=x[i][j]/contador;

for (i=1; i<=N; i++)
            for (j=1;j<=N; j++)
                fcorrelacion=fcorrelacion+x[i][j];

fcorrelacion=fcorrelacion/pow(N, 2);


// Escribimos los datos en un fichero.
fprintf(f1, "%e\t%e\n", temp, magnetizacion);
fprintf(f2, "%e\t%e\n", temp, energiamedia);
fprintf(f3, "%e\t%e\n", temp, caloresp);
fprintf(f4, "%e\t%e\n", temp, fcorrelacion);


}

printf("Proceso finalizado.");

fclose(f1);
fclose(f2);
fclose(f3);
fclose(f4);

return 0;

}
