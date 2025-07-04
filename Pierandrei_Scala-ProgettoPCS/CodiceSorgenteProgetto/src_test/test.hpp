 #pragma once

#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include "Utils.hpp"
#include "Poliedro.hpp"
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;
namespace LibreriaPoliedro {

TEST(CostruttoreTest, AreaFacce) {
	Poliedro T = Poliedro(3);
	Poliedro O = Poliedro(4);
	Poliedro I = Poliedro(5);
	
	bool t = true;
	bool o = true;
	bool i = true;
	
	double S = 0.0;
	
	for (int i = 0; i < T.NumFaccePoliedro; i++) {
		Vector3d A = T.CoordinateVertici.col(T.VerticiFacce[i][0]);
		Vector3d B = T.CoordinateVertici.col(T.VerticiFacce[i][1]);
		Vector3d C = T.CoordinateVertici.col(T.VerticiFacce[i][2]);
		S = 0.5 * ((B - A).cross(C - A)).norm();
		if (S == 0.0) {
			t = false;
		}
	}
	
	for (int i = 0; i < O.NumFaccePoliedro; i++) {
		Vector3d A = O.CoordinateVertici.col(O.VerticiFacce[i][0]);
		Vector3d B = O.CoordinateVertici.col(O.VerticiFacce[i][1]);
		Vector3d C = O.CoordinateVertici.col(O.VerticiFacce[i][2]);
		S = 0.5 * ((B - A).cross(C - A)).norm();
		if (S == 0.0) {
			o = false;
		}
	}
	
	for (int i = 0; i < O.NumFaccePoliedro; i++) {
		Vector3d A = I.CoordinateVertici.col(I.VerticiFacce[i][0]);
		Vector3d B = I.CoordinateVertici.col(I.VerticiFacce[i][1]);
		Vector3d C = I.CoordinateVertici.col(I.VerticiFacce[i][2]);
		S = 0.5 * ((B - A).cross(C - A)).norm();
		if (S == 0.0) {
			i = false;
		}
	}
	
	EXPECT_TRUE(t);
	EXPECT_TRUE(o);
	EXPECT_TRUE(i);
}
	

TEST(TriangolazioneTest, PoliedroTriangolato) {	
	// inizializzazione di un tetraedro, un ottaedro e un icosaedro
	Poliedro T = Poliedro(3);
	Poliedro O = Poliedro(4);
	Poliedro I = Poliedro(5);
	
	// triangolazione dei 3 poliedri con diversi passi
	T = Triangolato(T, 5);
	O = Triangolato(O, 4);
	I = Triangolato(I, 3);
	
	// verifica della corretta implementazione della funzione attraverso i numeri di vertici, spigoli e facce ATTESI	
	EXPECT_EQ(T.NumVerticiPoliedro, 52);
	EXPECT_EQ(O.NumVerticiPoliedro, 66);
	EXPECT_EQ(I.NumVerticiPoliedro, 92);
	
	EXPECT_EQ(T.NumSpigoliPoliedro, 150);
	EXPECT_EQ(O.NumSpigoliPoliedro, 192);
	EXPECT_EQ(I.NumSpigoliPoliedro, 270);
	
	EXPECT_EQ(T.NumFaccePoliedro, 100);
	EXPECT_EQ(O.NumFaccePoliedro, 128);
	EXPECT_EQ(I.NumFaccePoliedro, 180);
	

}

TEST(AdiacenzaTest, GrafoSemplice) {
    // inizializzazione di un quadrilatero come struttura poliedro
    Poliedro P;
    P.IdVertici = {0, 1, 2, 3};
    P.IdSpigoli = {0, 1, 2, 3};

    P.CoordinateVertici = MatrixXd(3, 4);
    P.CoordinateVertici <<
                -3, -3, 0, 3,
                -2, 2, 2, -2,
                0, 0, 0, 0;

    P.EstremiSpigoli = MatrixXi(2, 4);
    P.EstremiSpigoli <<
				0, 1, 2, 3,
				1, 2, 3, 0;				
				
	vector<vector<pair<unsigned int, double>>> A = Adiacenza(P);
	
	// verifica della corretta creazione del grafo di adiacenza 
	ASSERT_EQ(A.size(), 4);
	
	EXPECT_EQ(A[0][0].first, 1);
	EXPECT_DOUBLE_EQ(A[0][0].second, 4.0);
	EXPECT_EQ(A[0][1].first, 3);
	EXPECT_DOUBLE_EQ(A[0][1].second, 6.0);
	
	EXPECT_EQ(A[1][0].first, 0);
	EXPECT_DOUBLE_EQ(A[1][0].second, 4.0);
	EXPECT_EQ(A[1][1].first, 2);
	EXPECT_DOUBLE_EQ(A[1][1].second, 3.0);
	
	EXPECT_EQ(A[2][0].first, 1);
	EXPECT_DOUBLE_EQ(A[2][0].second, 3.0);
	EXPECT_EQ(A[2][1].first, 3);
	EXPECT_DOUBLE_EQ(A[2][1].second, 5.0);
	
	EXPECT_EQ(A[3][0].first, 2);
    EXPECT_DOUBLE_EQ(A[3][0].second, 5.0);	
	EXPECT_EQ(A[3][1].first, 0);
	EXPECT_DOUBLE_EQ(A[3][1].second, 6.0);
	
}
	
	
TEST(PercorsoMinimoTest, PercorsoMinimoSemplice) {
    
	// inizializzazione di un quadrilatero come struttura poliedro
	Poliedro P;
    P.IdVertici = {0, 1, 2, 3};
    P.IdSpigoli = {0, 1, 2, 3};

    P.CoordinateVertici = Eigen::MatrixXd(3, 4);
    P.CoordinateVertici <<
                -3, -3, 0, 3,
                -2, 2, 2, -2,
                0, 0, 0, 0;

    P.EstremiSpigoli = MatrixXi(2, 4);
    P.EstremiSpigoli <<
				0, 1, 2, 3,
				1, 2, 3, 0;

    std::vector<double> buffer;
    unsigned int origine = 0, estremo = 2;

    PercorsoMinimo(P, origine, estremo, buffer);
	
	// verifica del corretto aggiornamento del buffer (contenente le informazioni relative al percorso minimo trovato)
    ASSERT_EQ(buffer.size(), 4);
    EXPECT_EQ(buffer[0], 2.0);
    EXPECT_EQ(buffer[1], 2.0);
	EXPECT_EQ(buffer[2], 0.0);
    EXPECT_EQ(buffer[3], 0.0);
}

}