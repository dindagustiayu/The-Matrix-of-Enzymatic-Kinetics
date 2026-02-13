{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "da987ef1-c82a-4fdc-a09a-2c57828239f1",
   "metadata": {},
   "source": [
    " ---\n",
    "title: \"Michael-Menten Kinetics and Linear-Burk Plot\"\n",
    "date: \"2026-2-12\"\n",
    "categories: [Python 3, Jupyter Notebook, Numpy, Matplotlib, Enzymatic Kinetics]\n",
    "---"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "0011910b-f827-4a0a-8ba6-d1908bbbb675",
   "metadata": {},
   "source": [
    "[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "1f1976c2-22de-4136-bcb0-1e771f1993c5",
   "metadata": {},
   "source": [
    "# The Michael-Menten model and Linear-Burk Plot\n",
    "\n",
    "## Michael-Mentel Kinetics\n",
    "Michaelis-Menten kinetics describes the rate of enzymatic reactions, providing a model to understand how enzymes interact with substrates. The general reaction scheme of an enzyme-catalyzed reaction is as follows:\n",
    "<p align='center'>\n",
    "    $$E + S\\; \\overset{k_1}{\\underset{k_{-1}}{\\longrightleftharpoons}} \\;ES\\; \\overset{k_2}{\\longrightarrow}\\;P + E$$\n",
    "</p>\n",
    "The enzyme interacts with the substrate by binding to its active site to form the enzyme-substrate complex, ES. That reastion is followed by the decomposition of ES to generate the free enzyme, E, and the new product, P.\n",
    "\n",
    "| Rate constant | Reaction |\n",
    "| ------------- | -------- |\n",
    "| $k_{1}$ | The binding of the enzyme to the substrate forming the enzyme substrate complex. |\n",
    "| $k_{2}$ | Catalytic rate; the catalysis reaction producing the final reaction product and regenerating the free enzyme. This is the rate limiting step. |\n",
    "\n",
    "## Lineweaver-Burk Plot\n",
    "For many such catalyzed reactions, it is found that:\n",
    "<p align='center'>\n",
    "    $$v_{0}=\\frac{(V_{max}[S])}{(K_{M}+[S])}$$\n",
    "    $$\\frac{1}{v}=\\frac{(K_M + [S])}{v_{max}[S]}$$\n",
    "    $$\\frac{1}{v}=\\left(\\frac{K_M}{v_{max}} \\right)\\;\\left(\\frac{1}{[S]} \\right) + \\frac{1}{v_{max}}$$\n",
    "</p>\n",
    "Apply this to equation for a straight line $y = mx + b$ and we have:\n",
    "<p align='center'>\n",
    "    $$y=\\frac{1}{v}$$\n",
    "    $$x=\\frac{1}{[S]}$$\n",
    "    $$m=slope=\\frac{K_M}{V_{max}}$$\n",
    "    $$b=y-intercept=\\frac{1}{V_{max}}$$\n",
    "</p>\n",
    "When we plot $y=\\frac{1}{v}$ versus $x=\\frac{1}{[S]}$, we obtain a straight line.\n",
    "<p align='center'>\n",
    "    $$x-intercept=\\frac{-1}{K_M}$$\n",
    "    $$y-intercept=\\frac{1}{V_{max}}$$\n",
    "    $$slope=\\frac{K_M}{V_{max}}$$\n",
    "</p>\n",
    "\n",
    "![Figure. An example of a Lineweaver-Burke plot from Wikipedia](https://chem.libretexts.org/@api/deki/files/54231/2000px-Lineweaver-Burke_plot.svg.png?revision=1)\n",
    "\n",
    "## P19.5 Exercise \n",
    "Pepsin is one of the  main enzymes in the digestive systems of mammals, where it helps to break down proteins into smallers peptides which can be absorbed by the small intestine. In a study of the enzymatic kinetics of pepsin with the protein bovine serum albumin (S), the rate of reaction, $v$, was measured as a function of substrate concentration, [S]; these data are given in the file [pepsin.txt]()\n",
    "\n",
    "```\n",
    "    # Reaction rates, v (in mM.s-1), for enzymatic kinetics of pepsin with\n",
    "    # bovine serum albumin (S, in mM)\n",
    "\n",
    "    S / mM\tv / mM.s-1\n",
    "    0.1\t\t0.00339\n",
    "    0.2\t\t0.00549\n",
    "    0.5\t\t0.00875\n",
    "    1.0\t\t0.01090\n",
    "    5.0\t\t0.01358\n",
    "    10.0\t0.01401\n",
    "    20.0\t0.01423\n",
    "```\n",
    "\n",
    "In the case of a straight-line regression (\"simple linear regression\") on n points using the model function $y=a+b_x$ as follow:\n",
    "<p align='left'>\n",
    "    $$\\hat{\\beta}=\\begin{pmatrix} a\\\\b\\end{pmatrix}\\;where\\;a=\\frac{S_y S_{xx} - S_{xy} S_x}{nS_{xx} - S^2_x}\\;and\\;b=\\frac{nS_{xy} S_y S_x}{nS_{xx} - S^2_x},$$\n",
    "</p>\n",
    "and the summary statistics are defined as:\n",
    "<p align='center'>\n",
    "    $$S_x=\\overset{n}{\\underset{i=1}{\\sum x_i}},\\;S_y=\\overset{n}{\\underset{i=1}{\\sum Y_i}},\\;S_{xx}=\\overset{n}{\\underset{i=1}{\\sum x_i^2}},\\;S_{xy}=\\overset{n}{\\underset{i=1}{\\sum x_iy_i}}.$$\n",
    "</p>\n",
    "\n",
    "### Step 1. Plotted as a function of reaction rates for enzymatic kinetics of pepsin"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 82,
   "id": "c386f270-b78b-468c-83d3-b0a5ac27a6ba",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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M+SjvSQPPvigzJgxQ+bNm2eCkFY+tKLiDBBly5aV/v37m1DiPDYdprn11lulRo0aUqVKFYmKinK9j4anrVu3FvnnAJATFREAhQ4ir776ap7PawVj4cKFpm/kzTffNEMdWhm5//77zXBHQWh1I7fhjp9++slMkV28eLEJD84m0bVr18revXulQoUK5rHmzZtLUlKS+frHH38079+wYUPz/L59+8ythpnt27ebnylzNUTVqlXLhBB9Da3qXA/thdH30mGuyMhIWbZsmTkOAFlREQFQYPrhfeLECenUqVO++2mPiPZcaGjRD2atTmgg0YBSENojolWM3CoiOmziDCHq2LFjptrhDCHOx3QIx1kR0aChDakaSHTo5qmnnpJHHnlEbrrpJvO8sz/EqVy5cub2zz//vO5j/vrrr+XXX38136PnihAC5I4gAqDANFjExcW5PqDzowFC99XZJ5s3bzYzbSZPnlyg99Ogcfbs2VwrIi1btswRTjL3rVy8eFEOHDhghl3Ujh07XEFDG3J11o4Ou+gxXbp0Sfbs2ZMjiDiba6tWrVqg4wZwbQQRAAW2fPlyU4koDG0m1eGKgmjatKkZaslMp9PqkI0+53T06FETGjI/psEiIyPDhI6ff/7ZNNE6h140nOhwjE7P1Zk72kuis3yyD81o1USHVzJXXgC4Bz0iAHJITU2VQ4cOue4fPnzYVBp01olWOHR2Sn6LjyltGr3vvvtkyJAhZghEF/fSD/2ZM2fKPffcU6CzrsM6OqNFqyLaE+KshmijbOZ+Ez3GG264IctaI7qfTh3W9//qq68kKCjItc7JoEGDpFevXlK5cmVzX3s69PWdwzhO3377rXTu3JnfFMADCCIActDAoE2fTjrV1fnBfeedd5qhD13lND86Y0b3mz17tvznP/8xlQadiaLNq5MmTSrQWY+JiTGzZj799FOzQJozYGhPR+bhIR120cpHZrqfc1hGg4aGkMDAQHNfbzNXOfT5zNUU59CONppqzwcA9wtw6Oo+AHCddEimdevWZo2O7Eu86we+rudRqD9GAQHmA18rFLnRxczGjRtnhknymjLsCToFWIeidKEyAO5HjwiAAtEQorNScqNTYLUSkn3F0vyMGDHCfM+13HXXXaYacvLkSSlOWjV5/fXXi/U9AX9CRQSAW2hA0Gm2Sq8to70Y10OvR+O8josuJpZ52i2Ako8gAgAArGFoBgAAWEMQAQAA1hBEAACANQQRAABgDUEEAABYQxABAADWEEQAAIA1BBEAAGANQQQAAFhDEAEAANYQRAAAgNjyf6oHFzqUW6v/AAAAAElFTkSuQmCC",
      "text/plain": [
       "<Figure size 600x400 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "# Open file data\n",
    "S, v = np.genfromtxt('pepsin.txt', unpack=True, skip_header=3)\n",
    "\n",
    "# Transform data\n",
    "x = 1 / S   # 1/[S]\n",
    "y = 1 / v   # 1/v\n",
    "\n",
    "plt.figure(figsize=(6, 4))\n",
    "plt.plot(x, y, '+', markersize=10)\n",
    "plt.xlabel('$\\\\mathrm{1/[S]} \\\\; (mM)^{-1}$')\n",
    "plt.ylabel('$\\\\mathrm{1/[v]} \\\\; s(mM)^{-1}$')\n",
    "plt.savefig('Plotted as a function of enzymatic kinetics.svg', bbox_inches='tight')\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "6671bf73-9d90-4f0c-9f45-2305ff926ed9",
   "metadata": {},
   "source": [
    "## Step 2. Fitting a Line of Best fit\n",
    "The experiment were carried out at $35 ^{o}C$ and a pH of 2 with a totl concentration of pepsin of $[E]_{0}=0.028\\;mM$. Use `np.linalg.lstsq` to fit the data and obtain values for $K_M,\\;v_{max}\\;and\\;k_2$.\n",
    "\n",
    "For simple linear regression on $n$ points a straight line model, $y=a+bx$, \n",
    "<p align='center'>\n",
    "    $$\\hat{\\beta}=\\begin{pmatrix} a\\\\ b \\end{pmatrix}=(X^T X)^{-1} X^Ty$$\n",
    "</p>\n",
    "\n",
    "defines the best fit coefficients $a$ and $b$ where the matrices $y$ and $X$\n",
    "<p align='center'>\n",
    "    $$y=\\begin{pmatrix} y_1\\\\y_2\\\\ \\vdots \\\\y_n\\end{pmatrix},\\; and\\; X=\\begin{pmatrix}1 & x_1\\\\1 & x_2\\\\ \\vdots \\\\1 & x_n\\end{pmatrix}$$\n",
    "</p>\n",
    "We have\n",
    "<p align='center'>\n",
    "    $$X^T X=\\begin{pmatrix} 1 & 1 & 1 & \\cdots & 1\\\\ x_1 & x_2 & x_3 \\cdots & x_n\\end{pmatrix} \\; \\begin{pmatrix} 1 & x_1\\\\1 & x_2\\\\1 & x_3\\\\ \\vdots & \\vdots\\\\ 1 & x_n\\end{pmatrix} = \\begin{pmatrix} n & \\sum{x_i}\\\\ \\sum{x_i} & \\sum{x_i^2}\\end{pmatrix}$$\n",
    "</p>\n",
    "Taking the inverse,\n",
    "<p align='center'>\n",
    "    $$(X^TX X)^{-1}= \\frac{1}{\\Delta} \\;\\begin{pmatrix} \\sum{x_i^2} & -\\sum{x_i}\\\\ -\\sum{x_i} & n\\end{pmatrix}$$\n",
    "</p>\n",
    "where the determinant,\n",
    "<p align='center'>\n",
    "    $$\\Delta=n\\sum{x_i^2}-(\\sum{x_i})^2 = nS_{xx} - S_{x}^2$$\n",
    "</p>\n",
    "Next,\n",
    "\n",
    "<p align='center'>\n",
    "    $$\\begin{align*}(X^TX)^{-1}X^T&=\\frac{1}{\\Delta} \\begin{pmatrix} \\sum{x_i^2} & -\\sum{x_i}\\\\ -\\sum{x_i} & n\\end{pmatrix} \\begin{pmatrix} 1 & 1 & 1 & \\cdots & 1\\\\ x_1 & x_2 & x_3 & \\cdots & x_n\\end{pmatrix} \\\\&=\\frac{1}{\\Delta} \\begin{pmatrix} \\sum{x_i^2}-x_1 \\sum{x_1} & \\sum{x_i^2}-x_2 \\sum{x_i} & \\cdots & \\sum{x_i^2}-x_n \\sum{xi}\\\\ -\\sum{x_i}+ nx_1 & -\\sum{x_i}+nx_2 & \\cdots & -\\sum{x_i}+nx_n \\end{pmatrix} \\begin{pmatrix} y_1 \\\\ y_2 \\\\ \\vdots \\\\ y_n\\end{pmatrix} \\end{align*}$$\n",
    "</p>\n",
    "\n",
    "Finally,\n",
    "<p align='center'>\n",
    "    $$ \\begin{align*}\\hat{\\beta}&= (X^T X)^{-1}X^Ty \\\\ &= \\frac{1}{\\Delta} \\begin{pmatrix} \\sum{x_i^2}-x_1 \\sum{x_1} & \\sum{x_i^2}-x_2 \\sum{x_i} & \\cdots & \\sum{x_i^2}-x_n \\sum{xi}\\\\ -\\sum{x_i}+ nx_1 & -\\sum{x_i}+nx_2 & \\cdots & -\\sum{x_i}+nx_n \\end{pmatrix} \\begin{pmatrix} y_1 \\\\ y_2 \\\\ \\vdots \\\\ y_n\\end{pmatrix}=\\begin{pmatrix} a\\\\b \\end{pmatrix} \\end{align*}$$\n",
    "</p>\n",
    "where,\n",
    "<p align='center'>\n",
    "    $$\\begin{align} a&=\\frac{1}{\\Delta}[y_1(\\sum{x_i^2}-x_1 \\sum{x_i})+y_2(\\sum{x_i^2}-x_2\\sum{x_i})+ \\cdots + y_n (\\sum{x_i^2}-x_n\\sum{x_i})] \\\\&=\\frac{1}{\\Delta}[(\\sum{y_i}(\\sum{x_i^2})-(\\sum{x_i y_i})(\\sum{x_i})] \\\\ &=\\frac{S_yS_{xx}-S_{xy}S_x}{nS_{xx}-S_x^2} \\end{align}$$\n",
    "</p>\n",
    "and\n",
    "<p align='center'>\n",
    "    $$\\begin{align} b&= \\frac{1}{\\Delta}[y_1(nx_1-\\sum{x_i})+y_2(nx_2-\\sum{x_i})+ \\cdots + y_n(nx_n - \\sum{x_i})] \\\\&=\\frac{1}{\\Delta}[n(\\sum{x_iy_i}-(\\sum{y_i})(\\sum{x_i})]\\\\ &=\\frac{nS_{xy}-S_x S_y}{nS_{xx}-S_x^2}\\end{align}$$\n",
    "</p>"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 83,
   "id": "c7e1aa90-8668-4c42-999e-ec562f085104",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "vmax = 0.0145 mM/s\n",
      "Km   = 0.3267 mM\n",
      "k2   = 0.5166 s^-1\n"
     ]
    },
    {
     "data": {
      "image/png": "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",
      "text/plain": [
       "<Figure size 600x400 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "# load data from txt\n",
    "S, v = np.genfromtxt('pepsin.txt', unpack=True, skip_header=3)\n",
    "\n",
    "# Linear fit\n",
    "x = 1 / S\n",
    "y = 1 / v\n",
    "\n",
    "# Design matrix\n",
    "A = np.vstack([x, np.ones_like(x)]). T\n",
    "\n",
    "# Least squares fit\n",
    "m, b = np.linalg.lstsq(A, y, rcond=None)[0]\n",
    "\n",
    "# Parameters Lineweaver-Burk Plot\n",
    "vmax = 1 /b\n",
    "kM = m * vmax\n",
    "E0 = 0.028\n",
    "k2 = vmax /E0\n",
    "\n",
    "print(f\"vmax = {vmax:.4f} mM/s\")\n",
    "print(f\"Km   = {Km:.4f} mM\")\n",
    "print(f\"k2   = {k2:.4f} s^-1\")\n",
    "\n",
    "#Plot Lineweaver-Burk Plot\n",
    "x_fit = np.linspace(min(x), max(x), 100)\n",
    "y_fit = m * x_fit + b\n",
    "\n",
    "plt.figure(figsize=(6, 4))\n",
    "plt.plot(x, y, '+', color='blue', markersize=8, label='Data')\n",
    "plt.plot(x_fit, y_fit, color='black', label='Fit')\n",
    "plt.xlabel('$\\\\mathrm{1/[S]} \\\\; (mM)^{-1}$')\n",
    "plt.ylabel('$\\\\mathrm{1/[v]} \\\\; s(mM)^{-1}$')\n",
    "plt.title('Lineweaver-Burk Plot for Pepsin')\n",
    "plt.legend()\n",
    "plt.savefig('Linearweaver-Burk Plot for pepsin.svg', bbox_inches='tight')\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "d42b7338-725f-4029-9b3f-d57d8a0362f8",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.14.2"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
