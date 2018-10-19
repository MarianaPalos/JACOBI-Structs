#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

struct Jacobi
{
   float A[10][10], b[10], xo[10], x[10], D[10][10], L[10][10], U[10][10], T[10][10], R[10][10] ,C[10];
   int n,i,j,k,l,dom;

   float restbr[10],resmult[10],res[10], multAx[10], Ax[10], aux, B;

void Valores() //1. Ingresa los vlaores de los elementos de A,b, y xo
{
  cout << "Ingresa la dimensión N de la matriz cuadrada A: "<< endl; //Pide la dimension tanto de la matriz, como de los vectores
  cin>>n;

  cout <<  "Ingresa la matriz A: " << endl; //Pide los elementos de la matriz A.
   for(i=0;i<n;i++){
       for (j=0;j<n;j++){
           cout << "Ingresa el elemento a" << i+1 << j+1 << " de la matriz A: " << endl;
           cin >> A[i][j];
       }
   }

cout << "Ingresa los elementos de b: " << endl; //Elementos del arreglo b
for(i=0;i<n;i++){
        cout << "Ingresa el elemento b" << i+1 << " de la matriz b: " << endl;
        cin >> b[i];
}

cout << "Ingresa los elementos de xo: " << endl; //Elementos del arreglo xo
for(i=0;i<n;i++){
        cout << "Ingresa el elemento x" << i+1 << " de la matriz xo: " << endl;
        cin >> xo[i];
      }
}

void Dominante () //Convierte la matriz en diagonal dominante
{
  float suma, suma2,aux;
  int r=0,m;

  //Combirtiendola en diagonal Dominante
  for(i=0;i<n;i++)
  {
    suma=0;
    for(j=0;j<n;j++)
    {
      if (i!=j)
      {
        suma= suma + abs(A[i][j]);
      }
    }

    m=0;
    r=i;
    if(suma>abs(A[i][i]))//El renglon no cumple así que se cambia orden
    {
      for(int l=0;l<n;l++)//Este for solo es para que se repitan las lineas la n cantidad de veces.
      {//Aquí se hacen todas las posibles combinaciones para ver cual digito es el que iria en la diagonal
        suma2=0;
        aux=0;
        for(int k=0;k<n;k++)
        {
          if(k!=m)
            suma2= suma2 + abs(A[i][k]);
        }

        if(abs(A[i][m])>suma2)
        {
          for(i=0;i<n;i++)
          {
             aux=A[i][r];
             A[i][r]=A[i][m];
             A[i][m]=aux;
          }
        }
        m++;

      }

    }

  }

}

void Revisar (int op) //Funcion que revisa si una matriz es dominante
{
  float suma;
  int contador=0;

  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      suma=0;
      if (i!=j)
      {
        suma = suma + abs(A[i][j]);
      }
    }

    if(suma<abs(A[i][i]))
    {
        contador ++;
    //cout<<"\nRenglon " << i+1 << ": DOMINANTE"<< endl;


    if(op==1)
    {
      cout<< "\nLa matriz NO es dominante." << endl;
      Dominante();
      cout<<"\nLa Matriz despues del cambio es:"<<endl;
      Imprimir(A);
    }
    else
    {
      dom=0;
    }
  }
}


if (contador==3)
{
  cout <<"\nLa Matriz es DOMINANTE." << endl;
  dom=1;
  if (op==2)
  {
    Imprimir(A);
  }
}

}


void Imprimir (float M[][10]) //Imprime matrices cuadradas
{
  for(i=0;i<n;i++){
       cout << "|";
       for (j=0;j<n;j++){
           cout << " " << M[i][j] << " ";
       }
       cout << "|" << endl;
   }

}

void Imprimir2 (float M[10]) //Imprime arreglos de [1][n]
{
   for(i=0;i<n;i++){
           cout << "| " << M[i]<< " |" <<endl;
     }
   }

void Suma (float M1[][10], float M2[][10], float suma[][10]) //3. Funcion suma de dos matrices
{
  for(i = 0; i < n; i++)
  {
    for(j = 0; j < n; j++)
    {
      suma[i][j] = M1[i][j] + M2[i][j];
    }
  }
}

void RestaVect (float M1[10], float M2[10], float resta[10]) //3. Funcion suma de dos matrices
{
  for(i = 0; i < n; i++)
  {
      resta[i]= M1[i] - M2[i];
    }
  }

void Mult (float M1[][10], float M2[][10], float mult[][10]) //4. funcion multiplicacion de dos matrices
{

for(int i=0;i<n;i++){
  for(int j=0;j<n;j++){
          mult[i][j]=0;
          for(int k=0;k<6;k++){
  mult[i][j]= mult[i][j]+ M1[i][k] * M2[k][j];
  }
  }
}
}

void MultVector (float M1[][10], float V[10], float multi[10]) //5. Mult matriz por vector
{
  for(int i=0;i<n;i++){
    multi[i]=0;
    for(int j=0;j<n;j++){
    multi[i]= multi[i]+ M1[i][j] * V[j];
    }
    }
}

void InversaD ()//2. llena los elementos de la diagonal de D con 0, y los dem'as con 1/aij
{
for(i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
    {
      D[i][j]=0;
      if(j==i)
      D[i][j]=1.00/(A[i][j]);
    }
  }
}

void MatrizL()   //Saca la matriz L con la parte triangular inferior de A
{
  for(i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
    {
      L[i][j]=0;
      if (i>j)
      L[i][j]=A[i][j];
    }
  }
}


void MatrizU() //Saca la matriz U con la parte triangular superior de A
{
  for(i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
    {
      U[i][j]=0;
      if (i<j)
      U[i][j]=A[i][j];
    }
  }
}


//Metodo de Jacobi, usando las operaciones anteriores
void iteraciones() {
  int contador=0;
  aux=0;
  B=100;

  while (abs(aux-B)>0.0001)
  {
    contador++;
    for (i=0;i<n;i++)
    {
      x[i]=0;
      MultVector(R,xo,resmult);
      RestaVect(b,resmult,restbr);
      MultVector(D,restbr,res);
      x[i]=res[i];
      xo[i]=x[i];
      //cout << x[i]<<endl;
    }

    B=0;
    aux=0;

    for (i=0;i<n;i++)
    {
      MultVector(A,x,multAx);
    }
    aux = aux + multAx[i];
    B = B + b[i];
  }

cout << "Tu vector solución es: "<< endl;
for (i=0; i<n; i++)
{
  cout << "| ";
  cout << x[i] << " ";
  cout << " |" << endl;
}

cout << "Se hicieron " << contador << " iteraciones." << endl;
}

};


int main (void)
{
  Jacobi jac;

  jac.Valores();

  cout << "La Matriz A es: " <<": " << endl;
  jac.Imprimir(jac.A);

  cout << "b es : " <<": " << endl;
  jac.Imprimir2(jac.b);

  jac.Revisar(1);

  if(jac.dom!=1)
  {
    jac.Revisar(2);
  }

if (jac.dom==1)
{
  jac.InversaD();
  cout << "La Matriz D es: " <<": " << endl;
  jac.Imprimir(jac.D);

  jac.MatrizL();
  cout << "La Matriz L es: " <<": " << endl;
  jac.Imprimir(jac.L);

  jac.MatrizU();
  cout << "La Matriz U es: " <<": " << endl;
  jac.Imprimir(jac.U);

  jac.MultVector(jac.D,jac.b,jac.C);
  cout << "El Vector C es: " <<": " << endl;
  jac.Imprimir2(jac.C);

  jac.Suma(jac.L,jac.U,jac.R);
  cout << "La Matriz R resultante de L+U es: " <<": " << endl;
  jac.Imprimir(jac.R);

  jac.Mult(jac.D,jac.R,jac.T);
  cout << "La Matriz T es: " <<": " << endl;
  jac.Imprimir(jac.T);

  //jac.iteraciones();

  int contador=0, i;
  jac.aux=0;
  jac.B=100;

  while (abs(jac.aux-jac.B)>0.0000001)
  {
    contador++;
    for (i=0;i<jac.n;i++)
    {
      jac.x[i]=0;
      jac.MultVector(jac.R,jac.xo,jac.resmult);
      jac.RestaVect(jac.b,jac.resmult,jac.restbr);
      jac.MultVector(jac.D,jac.restbr,jac.res);
      jac.x[i]=jac.res[i];
      jac.xo[i]=jac.x[i];
      cout << jac.x[i]<<endl;
    }


    jac.B=0;
    jac.aux=0;

    for (i=0;i<jac.n;i++)
    {

    jac.MultVector(jac.A,jac.x,jac.multAx);

    jac.aux = jac.aux + jac.multAx[i];
    jac.B = jac.B + jac.b[i];
  }
  }

cout << "Tu vector solución es: "<< endl;
for (i=0; i<jac.n; i++)
{
  cout << "| ";
  cout << jac.x[i] << " ";
  cout << " |" << endl;
}

cout << "Se hicieron " << contador << " iteraciones." << endl;

}

else
cout << "\nNo se puede aplicar el metodo. No tiene convergencia." <<endl;


  return 0;

}
