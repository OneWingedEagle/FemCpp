//============================================================================
// Name        : FemCpp.cpp
// Author      : Hassan
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Vect.h"
using namespace std;

ostream &operator<<(ostream &os,  const Vect &v)  {

	for(int i=0;i<v.getLength();i++)
	os << v.get(i) << "\t" ;
	os << endl;

return os;
}

istream &operator>>(istream &is,  Vect &v) {

	double a;
	for(int i=0;i<3;i++){
		 is>>a;
		 v.set(a,i);

	}

return is;
}


int main() {
	cout << "!!!YY Hello World!!!" << endl; // prints !!!Hello World!!!


	Vect v1(3,4,5);
	Vect v2(4,2,0);
	Vect v3=v1+v2;
	// v3.show();
	//Vect v4(0,0,0);
	// cin>>v3;

	 cout<<v3;

	 //v3.hshow();

	return 0;
}
