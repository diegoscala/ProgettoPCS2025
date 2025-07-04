#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <queue>
#include <limits>
#include <utility>
#include <algorithm>
#include "UCDUtilities.hpp"
#include "Poliedro.hpp"


using namespace std;
using namespace Gedim;

namespace LibreriaPoliedro 
{

void Controllo(unsigned int& p, unsigned int& q, unsigned int& b, unsigned int& c, bool& G);

void Proiezione(Poliedro& Poliedro_Proiettato); 

Poliedro Duale(Poliedro& Originale);

unsigned int aggiunta_vertice(MatrixXd& CoordinateVertici, vector<unsigned int>& IdVertici, unsigned int& contatore_idvertici, const Vector3d& coord);

void aggiunta_spigolo(MatrixXi& EstremiSpigoli, vector<unsigned int>& IdSpigoli, unsigned int& contatore_idspigoli, map<std::pair<unsigned int, unsigned int>, unsigned int>& spigoli_id_map, unsigned int v1, unsigned int v2);

Poliedro Triangolato(const Poliedro& input, const unsigned int& b); 

void PolopolaSpigoliFacce(Poliedro& poliedro);

void EsportaCelle(const Poliedro& poliedro);

vector<vector<pair<unsigned int, double>>> Adiacenza(const Poliedro& poliedro);

void PercorsoMinimo(Poliedro& P, unsigned int& origine, unsigned int& estremo, vector<double>& bufferSpigoli);

void Colorazione(const std::vector<double>& bufferSpigoli, vector<UCDProperty<double>>& ColoreSpigoli);

}