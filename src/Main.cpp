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
Vect operator-(Vect  a){
	return a.times(-1);
}

bool operator==(Vect  a, Vect b){
	return a[0]==b[0];
}

int main() {
	cout << "!!!YY Hello World!!!" << endl; // prints !!!Hello World!!!


	Vect v1(3,4,5),	 v2(3,2,-5);
	//-v1;
	 ++v1;
	 v1.show();
	/*bool b=v1==v2;
	cout<<b<<endl;*/

	return 0;
}
