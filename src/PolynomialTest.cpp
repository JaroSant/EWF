//=======================================================================
// Copyright (C) 2003-2013 William Hallahan
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge,
// publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//=======================================================================

// PolynomialTest.cpp : Defines the entry point for the console application.
//

#include "Polynomial.h"
#include <iostream>
#include <vector>

void DisplayPolynomial(const Polynomial &polynomial);
void GetPolynomial(Polynomial &polynomial, int polynomial_type);

//======================================================================
//  Start of main program.
//======================================================================

int main(int argc, char *argv[]) {
  //------------------------------------------------------------------
  //  Get the type of test.
  //------------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "1. Find roots of the polynomial." << std::endl;
  std::cout << "2. Evaluate the polynomial at a real value" << std::endl;
  std::cout << "3. Evaluate the polynomial and its derivative at a real value"
            << std::endl;
  std::cout << "4. Evaluate the polynomial at a complex value" << std::endl;
  std::cout
      << "5. Evaluate the polynomial and its derivative at a complex value"
      << std::endl;
  std::cout << "6. Test polynomial arithmetic." << std::endl;
  std::cout << "7. Test polynomial division." << std::endl;
  std::cout << std::endl;
  std::cout << "Enter the type of test > ";
  int test_type;
  std::cin >> test_type;

  //------------------------------------------------------------------
  //  Get the type of polynomial.
  //------------------------------------------------------------------

  std::cout << "1. Arbitrary polynomial" << std::endl;
  std::cout
      << "2. Polynomial with maximum power and scalar value 1.0, the rest 0.0."
      << std::endl;
  std::cout << "3. Polynomial with  all coefficient equal to 1.0." << std::endl;
  std::cout << std::endl;
  std::cout << "Enter the type of polynomial > ";
  int polynomial_type;
  std::cin >> polynomial_type;
  std::cout << std::endl;

  //------------------------------------------------------------------
  //  Get a polynomial.
  //------------------------------------------------------------------

  Polynomial polynomial;

  GetPolynomial(polynomial, polynomial_type);

  //------------------------------------------------------------------
  //  Perform different processing for the different tests.
  //------------------------------------------------------------------

  switch (test_type) {
  case 1: {
    //----------------------------------------------------------
    //  Find the roots of the polynomial.
    //----------------------------------------------------------

    std::vector<double> real_vector;
    std::vector<double> imag_vector;

    int degree = polynomial.Degree();

    real_vector.resize(degree);
    imag_vector.resize(degree);

    double *real_vector_ptr = &real_vector[0];
    double *imag_vector_ptr = &imag_vector[0];

    int root_count = 0;

    if (polynomial.FindRoots(real_vector_ptr, imag_vector_ptr, &root_count) ==
        PolynomialRootFinder::SUCCESS) {
      int i = 0;

      for (i = 0; i < root_count; ++i) {
        std::cout << "Root " << i << " = " << real_vector_ptr[i] << " + i "
                  << imag_vector_ptr[i] << std::endl;
      }
    } else {
      std::cout << "Failed to find all roots." << std::endl;
    }
  }

  break;

  case 2: {
    //----------------------------------------------------------
    //  Evaluate the polynomial at a real value.
    //----------------------------------------------------------

    while (true) {
      std::cout << "Enter value > ";
      double xr;
      std::cin >> xr;
      std::cout << "P(" << xr << ") = " << polynomial.EvaluateReal(xr)
                << std::endl;
      std::cout << std::endl;
    }
  }

  break;

  case 3: {
    //----------------------------------------------------------
    //  Evaluate the polynomial and its derivative at a
    //  real value.
    //----------------------------------------------------------

    while (true) {
      std::cout << "Enter value > ";
      double xr;
      std::cin >> xr;

      double dr;
      double pr = polynomial.EvaluateReal(xr, dr);

      std::cout << "P(" << xr << ") = " << pr << std::endl;
      std::cout << "D(" << xr << ") = " << dr << std::endl;
      std::cout << std::endl;
    }
  }

  break;

  case 4: {
    //----------------------------------------------------------
    //  Evaluate the polynomial at a complex value.
    //----------------------------------------------------------

    while (true) {
      std::cout << "Enter real value > ";
      double xr;
      std::cin >> xr;

      std::cout << "Enter imaginary value > ";
      double xi;
      std::cin >> xi;

      double pr;
      double pi;

      polynomial.EvaluateComplex(xr, xi, pr, pi);

      std::cout << "P(" << xr << " + i " << xi << ") = " << pr << " + i " << pi
                << std::endl;
      std::cout << std::endl;
    }
  }

  break;

  case 5: {
    //----------------------------------------------------------
    //  Evaluate the polynomial and its derivative at a
    //  complex value.
    //----------------------------------------------------------

    while (true) {
      std::cout << "Enter real value > ";
      double xr;
      std::cin >> xr;

      std::cout << "Enter imaginary value > ";
      double xi;
      std::cin >> xi;

      double pr;
      double pi;
      double dr;
      double di;

      polynomial.EvaluateComplex(xr, xi, pr, pi, dr, di);

      std::cout << "P(" << xr << " + i " << xi << ") = " << pr << " + i " << pi
                << std::endl;
      std::cout << "D(" << xr << " + i " << xi << ") = " << dr << " + i " << di
                << std::endl;
      std::cout << std::endl;
    }
  }

  break;

  case 6: {
    //----------------------------------------------------------
    //  Test polynomial arithmetic.
    //  Test polynomial copy constructor and equals operator.
    //----------------------------------------------------------

    Polynomial p_0 = polynomial;
    Polynomial p_1;
    p_1 = p_0;

    //----------------------------------------------------------
    //  Test polynomial addition.
    //----------------------------------------------------------

    Polynomial p_sum = p_0 + p_1;

    std::cout << "The sum polynomial is:" << std::endl;
    std::cout << std::endl;
    DisplayPolynomial(p_sum);
    std::cout << std::endl;

    //----------------------------------------------------------
    //  Test polynomial subtraction.
    //----------------------------------------------------------

    std::cout << "The difference polynomial is:" << std::endl;
    Polynomial p_diff = p_0 - p_1;
    std::cout << std::endl;
    DisplayPolynomial(p_diff);
    std::cout << std::endl;

    //----------------------------------------------------------
    //  Test polynomial multiplication.
    //----------------------------------------------------------

    std::cout << "The product polynomial is:" << std::endl;
    Polynomial p_product = p_0 * p_1;
    std::cout << std::endl;
    DisplayPolynomial(p_product);
    std::cout << std::endl;
  }

  break;

  case 7: {
    //----------------------------------------------------------
    //  Get another polynomial that will be the divisor.
    //----------------------------------------------------------

    std::cout << "Enter the divisor polynomial." << std::endl;

    Polynomial divisor_polynomial;
    GetPolynomial(divisor_polynomial, 1);

    Polynomial quotient_polynomial;
    Polynomial remainder_polynomial;

    polynomial.Divide(divisor_polynomial, quotient_polynomial,
                      remainder_polynomial);

    //----------------------------------------------------------
    //  Display the quotient polynomial.
    //----------------------------------------------------------

    std::cout << "The quotient polynomial is:" << std::endl;
    std::cout << std::endl;
    DisplayPolynomial(quotient_polynomial);
    std::cout << std::endl;

    //----------------------------------------------------------
    //  Display the remainder polynomial.
    //----------------------------------------------------------

    std::cout << "The remainder polynomial is:" << std::endl;
    std::cout << std::endl;
    DisplayPolynomial(remainder_polynomial);
    std::cout << std::endl;
  }

  break;

  default:

    std::cout << "Invalid test type" << std::endl;
    return -1;
    break;
  }

  return 0;
}

//======================================================================
//  Function to display a polynomial.
//======================================================================

void DisplayPolynomial(const Polynomial &polynomial) {
  int power = 0;

  for (power = polynomial.Degree(); power > 0; --power) {
    //--------------------------------------------------------------
    //  Display the coefficient if it is not equal to one.
    //--------------------------------------------------------------

    if (polynomial[power] != 1.0) {
      std::cout << polynomial[power];
    }

    //--------------------------------------------------------------
    //  If this is not the scalar value, then display the variable
    //  X.
    //--------------------------------------------------------------

    if (power > 0) {
      std::cout << " X";
    }

    //--------------------------------------------------------------
    //  If this is higher than the first power, then display the
    //  exponent.
    //--------------------------------------------------------------

    if (power > 1) {
      std::cout << "^" << power;
    }

    //--------------------------------------------------------------
    //  Add each term together.
    //--------------------------------------------------------------

    std::cout << " + ";
  }

  //------------------------------------------------------------------
  //  Display the polynomial's scalar value.
  //------------------------------------------------------------------

  std::cout << polynomial[power] << std::endl;

  return;
}

//======================================================================
//  Function: GetPolynomial
//======================================================================

void GetPolynomial(Polynomial &polynomial, int polynomial_type) {
  //------------------------------------------------------------------
  //  Get the polynomial degree.
  //------------------------------------------------------------------

  std::cout << "Enter the polynomial degree > ";
  int degree = 0;
  std::cin >> degree;
  std::cout << std::endl;

  //------------------------------------------------------------------
  //  Create a buffer to contain the polynomial coefficients.
  //------------------------------------------------------------------

  std::vector<double> coefficient_vector;

  coefficient_vector.resize(degree + 1);

  double *coefficient_vector_ptr = &coefficient_vector[0];

  //------------------------------------------------------------------
  //  Create the specified type of polynomial.
  //------------------------------------------------------------------

  int i = 0;

  switch (polynomial_type) {
  case 1:

    //--------------------------------------------------------------
    //  Create an arbitrary polynomial.
    //--------------------------------------------------------------

    for (i = 0; i <= degree; ++i) {
      std::cout << "coefficient[" << i << "] = ";
      double temp;
      ;
      std::cin >> temp;
      coefficient_vector_ptr[i] = temp;
    }

    std::cout << std::endl;
    break;

  case 2:

    //--------------------------------------------------------------
    //  Create a polynomial with the maximum degree and the scalar
    //  value coefficients equal to 1.0 and all other coefficients
    //  equal to zero.
    //--------------------------------------------------------------

    for (i = 1; i < degree; ++i) {
      coefficient_vector_ptr[i] = 0;
    }

    coefficient_vector_ptr[0] = 1.0;
    coefficient_vector_ptr[degree] = 1.0;

    break;

  case 3:

    //--------------------------------------------------------------
    //  Create a polynomial with all coefficients equal to 1.0.
    //--------------------------------------------------------------

    for (i = 0; i <= degree; ++i) {
      coefficient_vector_ptr[i] = 1.0;
    }

    break;

  default:

    std::cout << "Invalid polynomial type" << std::endl;
    exit(-1);
  }

  //------------------------------------------------------------------
  //  Create an instance of class Polynomial.
  //------------------------------------------------------------------

  polynomial.SetCoefficients(coefficient_vector_ptr, degree);

  return;
}
