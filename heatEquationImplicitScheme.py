{
  "nbformat": 4,
  "nbformat_minor": 0,
  "metadata": {
    "colab": {
      "name": "heatEquationImplicitScheme.ipynb",
      "provenance": []
    },
    "kernelspec": {
      "name": "python3",
      "display_name": "Python 3"
    }
  },
  "cells": [
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "3Afs2D2FdHG1"
      },
      "source": [
        "#Semillero de Computación Científica\r\n",
        "## Marzo 4 de 2021\r\n",
        "## Tema: Método Implícita de las diferencias Finitas aplicada a la Ecuación de Calor\r\n"
      ]
    },
    {
      "cell_type": "code",
      "metadata": {
        "id": "0EyDRhmHdvZ1"
      },
      "source": [
        "# Vamos a definir los parámetros del problema\r\n",
        "\r\n",
        "import numpy as np\r\n",
        "import matplotlib.pyplot as plt\r\n",
        "import matplotlib\r\n",
        "from numpy import linalg"
      ],
      "execution_count": 12,
      "outputs": []
    },
    {
      "cell_type": "code",
      "metadata": {
        "id": "lwcTIPayiYNU"
      },
      "source": [
        "# funcion para encontrar calor h=dx, k=dt\r\n",
        "def solu_eq_calor_implicit(alpha,t0,a,b,nt,nx,dt,f):\r\n",
        "  # alpha: constante de calor\r\n",
        "  # t0: tiempo inicial\r\n",
        "  # a: extremo izquierdo\r\n",
        "  # b: extremo derecho\r\n",
        "  # nt: número de muestras en t\r\n",
        "  # nx: número de muestras en x\r\n",
        "  # dt: tamaño de muestra en tiempo\r\n",
        "  # f: condicion inicial\r\n",
        "  L=b-a\r\n",
        "  dx=L/(nx-1)\r\n",
        "  gamma=alpha*dt/dx**2\r\n",
        "\r\n",
        "  tmax=nt*dt\r\n",
        "  t=np.arange(t0,tmax,dt)\r\n",
        "\r\n",
        "  #inicialice el vector w\r\n",
        "  w=np.zeros((nx,nt))\r\n",
        "  X = np.linspace(a,b,w.shape[0])\r\n",
        "\r\n",
        "  print(\"dimensiones de w\", w.shape)\r\n",
        "\r\n",
        "  # condiciones de frontera (Dirichlet homogéneas)\r\n",
        "  for j in range(nt):\r\n",
        "    w[0,j]=0\r\n",
        "    w[nx-1,j]=0\r\n",
        "\r\n",
        "  # vector de la derecha y solución interna\r\n",
        "  rhs = np.zeros(nx-1)\r\n",
        "\r\n",
        "  # elaboración de matriz A\r\n",
        "  # Matriz A inicialización y lleno\r\n",
        "  A= np.zeros([nx-1,nx-1])\r\n",
        "\r\n",
        "  # condición inicial\r\n",
        "  for i in range(nx):\r\n",
        "    w[i,0]=f(X[i])\r\n",
        "\r\n",
        "  # llenado de la matrz A\r\n",
        "  for i in range(0,nx-1,1):\r\n",
        "    A[i][i] = 1.0 + 2. * gamma\r\n",
        "    if i> 0 : A[i][i-1] = - gamma\r\n",
        "    if i<nx-2: A[i][i+1] = - gamma\r\n",
        "\r\n",
        "  print(\"gamma=\",gamma)\r\n",
        "  # for i in range(nx-1):\r\n",
        "  # print(\", A[%d,0]=%f6.1\"%(i, A[i,0]))\r\n",
        "\r\n",
        "  # ciclo sobre tiempos\r\n",
        "  for j in range(0,nt-1):\r\n",
        "    # vector del lado derecho\r\n",
        "    for i in range(nx-1):\r\n",
        "      rhs[i]=w[i+1,j]\r\n",
        "    rhs[0] += gamma*w[0,j]\r\n",
        "    rhs[nx-2] += gamma*w[nx-1,j]\r\n",
        "\r\n",
        "    # resuelva Aw_j+1=w_j\r\n",
        "    x = linalg.solve(A, rhs)\r\n",
        "\r\n",
        "    # rotación\r\n",
        "    for i in range(nx-2):\r\n",
        "      w[i+1, j + 1] = x[i]\r\n",
        "\r\n",
        "  return w"
      ],
      "execution_count": 18,
      "outputs": []
    },
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "se269T97kEmm"
      },
      "source": [
        "Prueba del algoritmo\r\n"
      ]
    },
    {
      "cell_type": "code",
      "metadata": {
        "id": "UzgBv6PmkFp2"
      },
      "source": [
        "alpha=1\r\n",
        "a=0\r\n",
        "b=1\r\n",
        "nt=40001\r\n",
        "nx=101\r\n",
        "dt=0.00001\r\n",
        "t0=0\r\n",
        "x0=0"
      ],
      "execution_count": 19,
      "outputs": []
    },
    {
      "cell_type": "code",
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "_Pzf4bd1kIBy",
        "outputId": "50867dea-2481-4e10-bf9c-d4412a4b5b7e"
      },
      "source": [
        "L=b-a\r\n",
        "# definición de la condición inicial\r\n",
        "def f(x):\r\n",
        "  return 6*np.sin(np.pi*x/L)\r\n",
        "  \r\n",
        "w=solu_eq_calor_implicit(alpha,t0,a,b,nt,nx,dt,f)"
      ],
      "execution_count": 20,
      "outputs": [
        {
          "output_type": "stream",
          "text": [
            "dimensiones de w (101, 40001)\n",
            "gamma= 0.1\n"
          ],
          "name": "stdout"
        }
      ]
    },
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "u_8O50GikkG1"
      },
      "source": [
        "Grafica de la Solución"
      ]
    },
    {
      "cell_type": "code",
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/",
          "height": 295
        },
        "id": "H5PDly-NkmGx",
        "outputId": "af595cde-1523-459f-b70f-ceed6aef724e"
      },
      "source": [
        "plt.xlabel(r\"$x$\")\r\n",
        "plt.ylabel(r'$calor$')\r\n",
        "plt.title(r'Curvas de calor. Solución numérica')\r\n",
        "\r\n",
        "X = np.linspace(a,b,w.shape[0])\r\n",
        "for j in range(0,nt,10000):\r\n",
        "  t=j*dt\r\n",
        "  plt.plot(X,w[:,j], label=str(round(t,1)) + ' seg')\r\n",
        "  plt.legend(labelspacing=1, title=\"tiempo\")\r\n",
        "plt.grid(True)\r\n",
        "plt.show()"
      ],
      "execution_count": 21,
      "outputs": [
        {
          "output_type": "display_data",
          "data": {
            "image/png": "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\n",
            "text/plain": [
              "<Figure size 432x288 with 1 Axes>"
            ]
          },
          "metadata": {
            "tags": [],
            "needs_background": "light"
          }
        }
      ]
    }
  ]
}