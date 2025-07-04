#pragma once
#include <vector>
#include <map>
#include <array>
#include <numeric>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace LibreriaPoliedro{


struct Poliedro{
	
	// Cell0Ds (vertici)
    vector<unsigned int> IdVertici;     // ID vertici
    MatrixXd CoordinateVertici;      // Coordinate (x,y,z) Vertici

    // Cell1Ds (spigoli)
    vector<unsigned int> IdSpigoli;      // ID Spigoli
    MatrixXi EstremiSpigoli;           // Estremi (IdOrigine , IdFine) spigoli

    // Cell2Ds (facce)
    vector<unsigned int> IdFacce;       // ID Facce
	vector<unsigned int> SpecificaFacce;
    vector< vector<unsigned int> > VerticiFacce;   // Per ogni faccia(elemento vettore) lista id vertici 
    vector< vector<unsigned int> > SpigoliFacce;	// Per ogni faccia(elemento vettore) lista id spigoli
	
	// Cell3Ds (Poliedro)
    vector<unsigned int> IdPoliedro;       // ID Poliedro
	unsigned int NumVerticiPoliedro;
	unsigned int NumSpigoliPoliedro;	
	unsigned int NumFaccePoliedro;
    vector< vector<unsigned int> > VerticiPoliedro;   // Per ogni poliedro(elemento vettore) lista id vertici 
    vector< vector<unsigned int> > SpigoliPoliedro;      // Per ogni poliedro(elemento vettore) lista id spigoli
	vector< vector<unsigned int> > FaccePoliedro;    // Per ogni poliedro(elemento vettore) lista id facce
	
	
	// COSTRUTTORE POLIEDRO DEFAULT
	Poliedro() = default;
	
	//COSTRUTTORE SOLIDI PLATONICI per p = 3
	Poliedro(unsigned int q) {
        if (q == 3) {  // Tetraedro
            IdVertici.resize(4);
            iota(IdVertici.begin(), IdVertici.end(), 0);
            CoordinateVertici = MatrixXd(3, 4);
            CoordinateVertici <<
                1, -1, -1, 1,
                1, -1, 1, -1,
                1, 1, -1, -1;

            IdSpigoli = {0, 1, 2, 3, 4, 5};
            EstremiSpigoli = MatrixXi(2, 6);
            EstremiSpigoli <<
				0, 0, 0, 1, 1, 2,
				1, 2, 3, 2, 3, 3;

            IdFacce = {0, 1, 2, 3};
            VerticiFacce = {
                {0, 1, 2},
                {0, 1, 3},
                {0, 2, 3},
                {1, 2, 3}
            };

            SpigoliFacce = {
                {0, 1, 3},
                {0, 2, 4},
                {1, 2, 5},
                {3, 4, 5}
            };

            IdPoliedro = {0};
            VerticiPoliedro = { IdVertici };
            SpigoliPoliedro = { IdSpigoli };
            FaccePoliedro = { IdFacce };
			
			SpecificaFacce.assign(4,3);
			NumVerticiPoliedro = 4;
			NumSpigoliPoliedro = 6;
			NumFaccePoliedro = 4;
			
			

        } else if (q == 4) {  // Ottaedro
            IdVertici.resize(6);
            iota(IdVertici.begin(), IdVertici.end(), 0);
            CoordinateVertici = MatrixXd(3, 6);
            CoordinateVertici <<
				1, -1, 0, 0, 0, 0,
				0, 0, 1, -1, 0, 0,
				0, 0, 0, 0, 1, -1;

            IdSpigoli = {0,1,2,3,4,5,6,7,8,9,10,11};
            EstremiSpigoli = MatrixXi(2, 12);
            EstremiSpigoli <<
				0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3,
				2, 3, 4, 5, 2, 3, 4, 5, 4, 5, 4, 5;

            IdFacce = {0,1,2,3,4,5,6,7};
            VerticiFacce = {
                {0,2,4}, {0,4,3}, {0,3,5}, {0,5,2},
                {1,4,2}, {1,3,4}, {1,5,3}, {1,2,5}
            };

            SpigoliFacce = {
                {0,2,8}, {1,2,10}, {1,3,11}, {0,3,9},
                {4,6,8}, {5,6,10}, {5,7,11}, {4,7,9}
            };

            IdPoliedro = {0};
            VerticiPoliedro = { IdVertici };
            SpigoliPoliedro = { IdSpigoli };
            FaccePoliedro = { IdFacce };
			
			SpecificaFacce.assign(8,3);
			NumVerticiPoliedro = 6;
			NumSpigoliPoliedro = 12;
			NumFaccePoliedro = 8;

        } else if (q == 5) {  // Icosaedro
            const double phi = (1.0 + sqrt(5.0)) / 2.0;
            IdVertici.resize(12);
            iota(IdVertici.begin(), IdVertici.end(), 0);
            CoordinateVertici = MatrixXd(3, 12);
            CoordinateVertici <<
                -1, 1, -1, 1, 0, 0, 0, 0, phi, phi, -phi, -phi,
				phi, phi, -phi, -phi, -1, 1, -1, 1, 0, 0, 0, 0,
				0, 0, 0, 0, phi, phi, -phi, -phi, -1, 1, -1, 1;

            IdSpigoli.resize(30);
            iota(IdSpigoli.begin(), IdSpigoli.end(), 0);
            EstremiSpigoli = MatrixXi(2, 30);
            EstremiSpigoli <<
				0, 0, 0,  0,  0, 1, 1, 1, 1, 2, 2, 2,  2,  2, 3, 3, 3, 3, 4, 4,  4, 5,  5, 6, 6,  6, 7,  7, 8, 10,
				1, 5, 7, 10, 11, 5, 7, 8, 9, 3, 4, 6, 10, 11, 4, 6, 8, 9, 5, 9, 11, 9, 11, 7, 8, 10, 8, 10, 9, 11;
				
            // Dati facce icosaedro
            IdFacce.resize(20);
            iota(IdFacce.begin(), IdFacce.end(), 0);
            VerticiFacce = {
                {0,11,5}, {0,5,1}, {0,1,7}, {0,7,10}, {0,10,11},
                {1,5,9}, {5,11,4}, {11,10,2}, {10,7,6}, {7,1,8},
                {3,4,9}, {3,2,4}, {3,6,2}, {3,8,6}, {3,9,8},
                {4,2,11}, {2,6,10}, {6,8,7}, {8,9,1}, {9,5,4}
            };

            SpigoliFacce = {
				{4, 22, 1},    
				{1, 5, 0},     
				{0, 6, 2},     
				{2, 27, 3},    
				{3, 29, 4},    
				{5, 21, 8},    
				{22, 20, 18},   
				{29, 12, 13},    
				{27, 23, 25},    
				{6, 7, 26},    
				{14, 19, 17},  
				{9, 10, 14},  
				{15, 11, 9}, 
				{16, 24, 15},  
				{17, 28, 16},   
				{10, 13, 20},  
				{11, 25, 12},   
				{24, 26, 23},  
				{28, 8, 7},   
				{21, 18, 19}   
			};

            IdPoliedro = {0};
            VerticiPoliedro = { IdVertici };
            SpigoliPoliedro = { IdSpigoli };
            FaccePoliedro = { IdFacce };
			
			SpecificaFacce.assign(20,3);
			NumVerticiPoliedro = 12;
			NumSpigoliPoliedro = 30;
			NumFaccePoliedro = 20;

        }
    } //fine costruttore

	
};


}
