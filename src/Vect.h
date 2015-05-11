/*
 * Vect.h
 *
 *  Created on: 2015/01/23
 *      Author: Hassan
 */

#ifndef VECT_H_
#define VECT_H_
using namespace std;
#include <iostream>
#include <vector>


class Vect{
public:
	Vect(){};

	Vect(double a, double b){
		length=2;
		el.resize(length);
		el[0]=a;
		el[1]=b;

	};

	Vect(double a, double b, double c){
		length=3;
		el.resize(length);
		el[0]=a;
		el[1]=b;
		el[2]=c;

	};

	Vect(double array[],int size){
		length=size;

		el.resize(length);

		for(int i=0;i<length;i++)
			el[i]=array[i];
	};

	Vect(int L){
		length=L;
		el.resize(length);

	};

	Vect add(Vect b){
		if(this->length!=b.length)
			cout<<"Vectors have different length";


		double c[this->length];
		for(int i=0;i<length;i++)
			c[i]=this->get(i)+b.get(i);


		return Vect(c,this->length);

	};


	Vect sub(Vect b){
		if(this->length!=b.length)
			cout<<"Vectors have different length";


		double c[this->length];
		for(int i=0;i<length;i++)
			c[i]=this->get(i)-b.get(i);


		return Vect(c,this->length);

	};

	Vect times(double a){

		double c[this->length];
		for(int i=0;i<length;i++)
			c[i]=this->get(i)*a;


		return Vect(c,this->length);

	};

	Vect copy(){

		double c[this->length];
		for(int i=0;i<length;i++)
			c[i]=this->get(i);


		return Vect(c,this->length);

	};


	Vect operator+(Vect  b){
		return this->add(b);
	}

	void operator+=(Vect  b){

		if(this->length!=b.length)
					cout<<"Vectors have different length";

		for(int i=0;i<length;i++)
		this->set(this->get(i)+b.get(i),i);
		}

	Vect operator-(Vect  b){
		return this->sub(b);
	}

	Vect operator*(double  a){
		return times(a);
	}

	void operator-(){
		for(int i=0;i<length;i++)
					el[i]=-el[i];
	}

	void operator++(){
		for(int i=0;i<length;i++)
					el[i]+=1;
	}

	double operator[](int k){
		return this->get(k);
	}



	Vect inv(){

		double c[this->length];
		for(int i=0;i<length;i++)
			c[i]=1.0/get(i);

		return Vect(c,this->length);
	}


	Vect aug(Vect  b){
		int L=this->length+b.length;
		double c[L];
		for(int i=0;i<L;i++)
			if(i<length)
				c[i]=get(i);
			else
				c[i]=b.get(i-length);
		return Vect(c,L);

	}


	int getLength(){
		return this->length;
	}


	double get(int i){
		return this->el[i];
	}


	void set(double a, int i){
		el[i]=a;
	}



	void hshow(){
		int L=this->length;
		for(int i=0;i<L;i++)
			cout<<get(i)<<"\t";
		cout<<endl;

	}

	void show(){
		int L=this->length;
		for(int i=0;i<L;i++)
			cout<<get(i)<<endl;

	}


private:
	int length;
	vector<double> el;

};



#endif /* VECT_H_ */
