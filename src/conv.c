     void convolve(double *x, int *nx, double *y, int *ny, double *xy)
     {
       int i, j, nxy = *nx + *ny - 1;
     
       for(i = 0; i < nxy; i++)
         xy[i] = 0.0;
       for(i = 0; i < *nx; i++)
         for(j = 0; j < *ny; j++)
           xy[i + j] += x[i] * y[j];
     }

/*
gcc -c -fPIC conv.c -o conv.o 
gcc -shared -Wl,-soname,libconv.so.1 -o libconv.so.1.0.1  conv.o
 */
