/**
 * @file weibull_nr.c
 * @author Francisco Ariel Vázquez Chávez (vazquez.francisco@colpos.mx)
 * @brief Programa para realizar el algoritmo NR para estimar parametros de la distribución weibull con censura por derecha tipo I y la matriz de varianza y covarianzas en C.
 * @version 0.1
 * @date 2022-12-20
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <rmath.h>
#include <R.h>


double nocens(int n, double *di)
{
    double s;
    s = 0;
    int i;
    for (i = 0; i < n; i++)
    {
        s = s + di[i];
    }
    return (s);
}

/// @brief Función auxiliar que calcula la suma de t_i^gama.
/// @param ti Double. Vector que contiene los tiempos de falla.
/// @param r Entero. Tamaño del vector.
/// @param gama Double. Parámetro de forma.
/// @return Double. Suma de t_i^gama.
double stg(double *ti, int r, double gama)
{
    double suma = 0;
    int i;
    for (i = 0; i < r; i++)
    {
        suma = suma + pow(ti[i], gama);
    }
    return (suma);
}

/// @brief Función auxiliar que calcula la suma de log(t_i).
/// @param r Entero. Tamaño del vector.
/// @param ti Double. Vector que contiene los tiempos de falla.
/// @param di Double. Vector que contiene la variable indicadora de censura.
/// @return Double. Suma de log(t_i).
double slt(int r, double *ti, double *di)
{
    double suma = 0;
    int i;
    for (i = 0; i < r; i++)
    {
        suma = suma + di[i] * log(ti[i]);
    }
    return (suma);
}

/// @brief Función auxiliar que calcula la suma de t_i^gama*log(t_i).
/// @param ti Double. Vector que contiene los tiempos de falla.
/// @param r Entero. Tamaño del vector.
/// @param gama Double. Parámetro de forma.
/// @return Double. Suma de t_i^gama*log(t_i).
double sltg(double *ti, int r, double gama)
{
    double suma = 0;
    int i;
    for (i = 0; i < r; i++)
    {
        suma = suma + pow(ti[i], gama) * log(ti[i]);
    }
    return (suma);
}

/// @brief Función auxiliar que calcula la suma de t_i^gama*log(t_i)^2.
/// @param ti Double. Vector que contiene los tiempos de falla.
/// @param r Entero. Tamaño del vector.
/// @param gama Double. Parámetro de forma.
/// @return Double. Suma de t_i^gama*log(t_i)^2.
double sltg2(double *ti, int r, double gama)
{
    double suma = 0;
    int i;
    for (i = 0; i < r; i++)
    {
        suma = suma + pow(ti[i], gama) * pow(log(ti[i]), 2);
    }
    return (suma);
}

/// @brief Calcula la log-verosimilitud para datos Weibull con censura por derecha tipo I.
/// @param beta Double. Parámetro de escala, debe ser mayor a cero.
/// @param gama Double. Parámetro de forma, debe ser mayor a cero.
/// @param n Entero. Dimensión del vector.
/// @param ti Double. Vector de dimensión n que contiene los tiempos de falla.
/// @param di Double. Vector de dimensión n que indica si la observación es completa (1) o censurada (0).
/// @return Double. Log-verosimilitud
double lvero(double beta, double gama, int n, double *ti, double *di)
{
    double lv;
    double r;
    r = nocens(n, di);
    lv = r * (log(gama) - gama * log(beta)) + (gama - 1) * slt(r, ti, di) - pow(1 / beta, gama) * stg(ti, n, gama);
    return (lv);
}

/// @brief Calcula el primer elemento del gradiente o Vector de primeras derivadas de la función de log-verosimilitud con respecto a beta.
/// @param beta Double. Parámetro de escala, debe ser mayor a cero.
/// @param gama Double. Parámetro de forma, debe ser mayor a cero.
/// @param n Entero. Tamaño de la muestra.
/// @param ti Double. Vector de dimensión n que contiene los tiempos de falla no censurados.
/// @param di Double. Vector de dimensión n que indica si la observación es completa (1) o censurada (0).
/// @return Double. Elemento S1 de la función Score
double s1(double beta, double gama, int n, double *ti, double *di)
{
    double s;
    double r;
    r = nocens(n, di);
    s = -1 * gama / beta * r + gama * pow(1 / beta, gama + 1) * stg(ti, n, gama);
    return s;
}

/// @brief Calcula el segundo elemento del gradiente o Vector de primeras derivadas de la función de log-verosimilitud.
/// @param beta Double. Parámetro de escala, debe ser mayor a cero.
/// @param gama Double. Parámetro de forma, debe ser mayor a cero.
/// @param n Entero. Tamaño de la muestra.
/// @param ti Double. Vector de dimensión n que contiene los tiempos de falla no censurados.
/// @param di Double. Vector de dimensión n que indica si la observación es completa (1) o censurada (0).
/// @return Double. Elemento S2 de la función Score
double s2(double beta, double gama, int n, double *ti, double *di)
{
    double s;
    double r;
    r = nocens(n, di);
    s = r * (1 / gama - log(beta)) + slt(n, ti, di) - pow(1 / beta, gama) * (sltg(ti, n, gama) + log(1 / beta) * stg(ti, n, gama));
    return s;
}

/// @brief Calcula el gradiente o Vector de primeras derivadas de la función de log-verosimilitud.
/// @param beta Double. Parámetro de escala, debe ser mayor a cero.
/// @param gama Double. Parámetro de forma, debe ser mayor a cero.
/// @param n Entero. Tamaño de la muestra.
/// @param ti Double. Vector de dimensión n que contiene los tiempos de falla no censurados.
/// @param di Double. Vector de dimensión n que indica si la observación es completa (1) o censurada (0).
/// @param u Double. Vector de dimensión 2 en donde se almacenará la información.
void score(double beta, double gama, int n, double *ti, double *di, double *u)
{
    u[0] = s1(beta, gama, n, ti, di);
    u[1] = s2(beta, gama, n, ti, di);
}

/// @brief Función auxiliar para calcular el elemento h11 de la matriz Hessiana.
/// @param beta Double. Parámetro de escala, debe ser mayor a cero.
/// @param gama Double. Parámetro de forma, debe ser mayor a cero.
/// @param n Entero. Tamaño de la muestra.
/// @param ti Double. Vector de dimensión n que contiene los tiempos de falla no censurados.
/// @param di Double. Vector de dimensión n que indica si la observación es completa (1) o censurada (0).
/// @return Double. Valor del elemento h11.
double h11(double beta, double gama, int n, double *ti, double *di)
{
    double r, s;
    r = nocens(n, di);
    s = gama / pow(beta, 2) * r - gama * (gama + 1) * pow(1 / beta, gama + 2) * stg(ti, n, gama);
    return (s);
}

/// @brief Función auxiliar para calcular el elemento h12 de la matriz Hessiana.
/// @param beta Double. Parámetro de escala, debe ser mayor a cero.
/// @param gama Double. Parámetro de forma, debe ser mayor a cero.
/// @param n Entero. Tamaño de la muestra.
/// @param ti Double. Vector de dimensión n que contiene los tiempos de falla no censurados.
/// @param di Double. Vector de dimensión n que indica si la observación es completa (1) o censurada (0).
/// @return Double. Valor del elemento h12.
double h12(double beta, double gama, int n, double *ti, double *di)
{
    double r, s;
    r = nocens(n, di);
    s = -1 * r / beta + pow(1 / beta, gama + 1) * (gama * sltg(ti, n, gama) + (1 + gama * log(1 / beta)) * stg(ti, n, gama));
    return (s);
}

/// @brief Función auxiliar para calcular el elemento h22 de la matriz Hessiana.
/// @param beta Double. Parámetro de escala, debe ser mayor a cero.
/// @param gama Double. Parámetro de forma, debe ser mayor a cero.
/// @param n Entero. Tamaño de la muestra.
/// @param ti Double. Vector de dimensión n que contiene los tiempos de falla no censurados.
/// @param di Double. Vector de dimensión n que indica si la observación es completa (1) o censurada (0).
/// @return Double. Valor del elemento h22.
double h22(double beta, double gama, int n, double *ti, double *di)
{
    double r, s;
    r = nocens(n, di);
    s = -1 * r / pow(gama, 2) - pow(1 / beta, gama) * (sltg2(ti, n, gama) + 2 * log(1 / beta) * sltg(ti, n, gama) + pow(log(1 / beta), 2) * stg(ti, n, gama));
    return (s);
}

/// @brief Calcula el hessiano o Matriz de segundas derivadas de la función de log-verosimilitud.
/// @param beta Double. Parámetro de escala, debe ser mayor a cero.
/// @param gama Double. Parámetro de forma, debe ser mayor a cero.
/// @param n Entero. Tamaño de la muestra.
/// @param ti Double. Vector de dimensión n que contiene los tiempos de falla no censurados.
/// @param di Double. Vector de dimensión n que indica si la observación es completa (1) o censurada (0).
/// @param h Double. Vector de dimensión 4 en donde se almacenará la información.
void hessiano(double beta, double gama, int n, double *ti, double *di, double *h)
{
    h[0] = h11(beta, gama, n, ti, di);
    h[1] = h12(beta, gama, n, ti, di);
    h[2] = h12(beta, gama, n, ti, di);
    h[3] = h22(beta, gama, n, ti, di);
}

/// @brief Calcula la matriz de covarianzas de la estimadores.
/// @param beta Double. Parámetro de escala, debe ser mayor a cero.
/// @param gama Double. Parámetro de forma, debe ser mayor a cero.
/// @param n Entero. Tamaño de la muestra.
/// @param ti Double. Vector de dimensión n que contiene los tiempos de falla no censurados.
/// @param di Double. Vector de dimensión n que indica si la observación es completa (1) o censurada (0).
/// @param v Double. Vector de dimensión 4 en donde se almacenará la información.
void VC(double beta, double gama, int n, double *ti, double *di, double *v)
{
    double det;
    double h_11, h_12, h_22;
    h_11 = h11(beta, gama, n, ti, di);
    h_12 = h12(beta, gama, n, ti, di);
    h_22 = h22(beta, gama, n, ti, di);
    det = h_11 * h_22 - pow(h_12, 2);
    v[0] = -1 * h_22 / (det);
    v[1] = h_12 / (det);
    v[2] = h_12 / (det);
    v[3] = -1 * h_11 / (det);
}

/// @brief Función Auxiliar que calcula Inv(H)*S.
/// @param x0 Double. Vector de dimensión 2 que contiene valores iniciales para los parámetros beta y gama.
/// @param n Entero. Tamaño de la muestra.
/// @param ti Double. Vector de dimensión n que contiene los tiempos de falla.
/// @param di Double. Vector de dimensión n que indica si la observación es completa (1) o censurada (0).
/// @param xfin Double. Vector de dimensión 2 en donde se almacenará la información.
void aux_nr(double *x0, int *n, double *ti, double *di, double *xfin)
{
    double beta, gama;
    beta = x0[0];
    gama = x0[1];
    int ene;
    ene = *n;
    double det, s_1, s_2;
    double h_11, h_12, h_22;
    h_11 = h11(beta, gama, ene, ti, di);
    h_12 = h12(beta, gama, ene, ti, di);
    h_22 = h22(beta, gama, ene, ti, di);
    det = h_11 * h_22 - pow(h_12, 2);
    s_1 = s1(beta, gama, ene, ti, di);
    s_2 = s2(beta, gama, ene, ti, di);
    xfin[0] = (h_22 * s_1 - h_12 * s_2) / det;
    xfin[1] = (h_11 * s_2 - h_12 * s_1) / det;
}
/// @brief Estima los parámetros de escala y forma para la distribución Weibull con censura por derecha tipo I usando el algoritmo de Newton-Raphson.
/// @param x0 Double. Vector de dimensión 2 que contiene valores iniciales para los parámetros beta y gama.
/// @param ops Double. Vector de dimensión 2 que contiene el número de iteraciones y la tolerancia.
/// @param n Entero. Tamaño de la muestra.
/// @param ti Double. Vector de dimensión n que contiene los tiempos de falla.
/// @param di Double. Vector de dimensión n que indica si la observación es completa (1) o censurada (0).
/// @param cod Entero. Regresa el código del algoritmo.
/// @param xfin Double. Vector de dimensión 2 en donde se almacenará la información.
/// @param vcov Double. Matriz en donde se guardará la matriz de covarianzas.
/// @param hist Entero. Indica si se debe imprimir el historial de iteraciones (1).
void weibull_nr(double *x0, double *ops, int *n, double *ti, double *di, int *cod, double *xfin, double *vcov, int *hist)
{
    int i = 1, codigo;
    double niter, tol;
    niter = ops[0];
    tol = ops[1];
    double xinc[2], xini[2];
    xini[0] = x0[0];
    xini[1] = x0[1];
    if (*hist == 1)
    {
        Rprintf("Valores iniciales: beta: %lf gama: %lf logv.: %lf\n", xini[0], xini[1], lvero(xini[0], xini[1], *n, ti, di));
    }

    while (i <= niter)
    {
        aux_nr(xini, n, ti, di, xinc);
        xfin[0] = xini[0] - xinc[0];
        xfin[1] = xini[1] - xinc[1];
        if (*hist == 1)
        {
            Rprintf("Paso %d: beta: %lf gama: %lf logv.: %lf\n", i, xfin[0], xfin[1], lvero(xfin[0], xfin[1], *n, ti, di));
        }
        if (fabs(xfin[0] - xini[0]) < tol && fabs(xfin[1] - xini[1]) < tol)
        {
            printf("El algoritmo ha convergido en %d pasos.\n", i);
            codigo = 1;
            *cod = codigo;
            VC(xfin[0], xfin[1], *n, ti, di, vcov);
        }
        xini[0] = xfin[0];
        xini[1] = xfin[1];
        i = i + 1;
        if (*cod == 1)
            break;
    }
    if (i >= niter && !(fabs(xfin[0] - xini[0]) < tol && fabs(xfin[1] - xini[1]) < tol))
    {
        printf("El algoritmo no encontró un óptimo.\n");
        codigo = -1;
        *cod = codigo;
    }
}

