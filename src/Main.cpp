//============================================================================
// Name        : Main.cpp
// Author      : Hassan
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Vect.h"
using namespace std;

int main() {
	cout << "!!!YY Hello World!!!" << endl; // prints !!!Hello World!!!


	Vect v1(3,4,5);
	Vect v2(4,2,0);
	 Vect v3=v1+v2;
	 v3.show();

	return 0;
}
