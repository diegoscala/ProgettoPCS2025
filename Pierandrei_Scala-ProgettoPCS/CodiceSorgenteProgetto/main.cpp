#include <iostream>
#include <sstream>
#include <string>
#include <queue>
#include <limits>
#include <utility>
#include <algorithm>
#include <vector>
#include "Poliedro.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"


using namespace std;
using namespace Eigen;
using namespace LibreriaPoliedro;

int main() {
    string input;
    unsigned int p = 0, q = 0, b = 0, c = 0, o = 0, e = 0; // {p,q} simbolo di schlafli, (b,c) parametri di triangolarizzazione, (o,e) eventuali indici di percorso minimo
    bool G = false; //G controlla l'uso del duale per generare un eventuale poliedro di Goldberg
    bool percorso_minimo = false; // controlla il caso in cui si debba calcolare il percorso minimo

	//ricezione dei dati in input
    cout << "Inserisci 4 o 6 numeri separati da spazio: ";
    getline(std::cin, input);
    istringstream iss(input);

    // lettura dei numeri 
    vector<unsigned int> nums;
    unsigned int temp;
    while (iss >> temp) {
        nums.push_back(temp);
    }
	// controllo quantità di input passati
    if (nums.size() != 4 && nums.size() != 6) {
        std::cerr << "Errore: devi inserire esattamente 4 o 6 numeri.\n";
        return 1;
    }

    // assegnazione dei valori
    p = nums[0];
    q = nums[1];
    b = nums[2];
    c = nums[3];
    if (nums.size() == 6) {
        o = nums[4];
        e = nums[5];
        percorso_minimo = true;
    }

    Controllo(p, q, b, c, G); // controllo validità e aggiornamento dei valori

    Poliedro P(q); // generazione solido platonico relativo a {p,q}
	
    P = Triangolato(P, b); //triangolarizzazione del solido platonico con passo b

	// se il solido platonico di partenza era un cubo o un dodecaedro occorre calcolare il poliedro di Goldberg
    if (G) {
		P = Duale(P); 
    }

    Proiezione(P); // proiettiamo i vertici del poliedro ottenuto sulla sfera unitario 
	
	// se sono stati inseriti anche gli id di inizio e fine percorso si calcola il percorso minimo fra questi due vertici e si esportano gli output
    if (percorso_minimo) {
        if (o < P.IdVertici.size() && e < P.IdVertici.size()) { // se gli id sono validi si procede col percorso minimo prima di esportare gli output
			
			// inizializzazione delle strutture per il colore degli spigoli che fanno parte del percorso minimo
            vector<UCDProperty<double>> ColoreSpigoli;
            vector<double> bufferSpigoli;
			 
			PercorsoMinimo(P, o, e, bufferSpigoli); // calcolo percorso minimo
			
			Colorazione(bufferSpigoli, ColoreSpigoli); // completamento della struttura UCDProperty

			// esportazione output per visualizzazione in paraview
            Gedim::UCDUtilities utilities;
            utilities.ExportPoints("./Vertici.inp", P.CoordinateVertici);
            utilities.ExportSegments("./Spigoli.inp", P.CoordinateVertici, P.EstremiSpigoli, {}, ColoreSpigoli);
        }
		else { // se gli id non sono validi stampa errore
            std::cerr << "Errore: gli ID dei vertici inseriti per il percorso minimo non sono validi.\n";
            return 1;
        }
    }
	else {
		// se non vengono passati gli id si esportano gli output per la visualizzazione in paraview
		Gedim::UCDUtilities utilities;
		utilities.ExportPoints("./Vertici.inp", P.CoordinateVertici);
		utilities.ExportSegments("./Spigoli.inp", P.CoordinateVertici, P.EstremiSpigoli);
	}
	
	EsportaCelle(P);

    return 0;
}
