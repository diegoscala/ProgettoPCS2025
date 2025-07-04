#include "Utils.hpp"

namespace LibreriaPoliedro {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// funzione che verifica la correttezza del simbolo di schlafli, dei parametri della triangolarizzazione e li processa
void Controllo(unsigned int& p, unsigned int& q, unsigned int& b, unsigned int& c, bool& G) {
	bool bc = ((b == 0 && c > 0) || (c == 0 && b > 0));
    bool presenza_tre = (p == 3 || q == 3);
    bool altro = ((p == 3 && (q == 3 || q == 4 || q == 5)) ||
                  (q == 3 && (p == 3 || p == 4 || p == 5)));
    bool pq = presenza_tre && altro;

    if (!bc || !pq) {
        std::cerr << "Errore nei valori:\n";
        if (!bc) {
            std::cerr << "- b e c devono essere uno uguale a 0 e l'altro maggiore di 0 (b = "
                      << b << ", c = " << c << ")\n";
        }
        if (!pq) {
            std::cerr << "- tra p e q, uno deve essere 3 e l'altro deve essere 3, 4 o 5 (p = "
                      << p << ", q = " << q << ")\n";
        }
    }

    // portiamo il valore 3 in p per chiamare il costruttore
    if (p != 3 && q == 3) {
        std::swap(p, q);
        G = true;
    }
	// portiamo il valore 0 in c per chiamare la triangolarizzazione
    if (c != 0 && b == 0) {
        std::swap(b, c);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// funzione per proiettare i vertici di un poliedro sulla sfera di raggio unitario
void Proiezione(Poliedro& Poliedro_Proiettato) {
    for (int i = 0; i < Poliedro_Proiettato.CoordinateVertici.cols(); ++i) {
        Poliedro_Proiettato.CoordinateVertici.col(i).normalize();
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
// funzione per calcolare il duale di un polidro
Poliedro Duale(Poliedro& Originale) {

	Poliedro Duale;


	//VERTICI
	unsigned int numFacce = Originale.IdFacce.size();
	Duale.IdVertici.resize(numFacce);
	Duale.CoordinateVertici = MatrixXd(3, numFacce);

	//Calcolo Baricentri
	for (size_t i=0 ; i < numFacce ; ++i){
		Vector3d baricentro = Vector3d::Zero();
		for (auto id_vertice : Originale.VerticiFacce[i]){
			baricentro += Originale.CoordinateVertici.col(id_vertice);
		}
		baricentro /= Originale.VerticiFacce[i].size();
		Duale.CoordinateVertici.col(i) = baricentro;
	}

	// Nuovi vertici (ID semplici da 0 a numFacce-1)
    for (size_t i = 0; i < Originale.IdFacce.size(); ++i)
		Duale.IdVertici[i] = i;

	//SPIGOLI

	// Mappa che associa a ciascuno spigolo le facce che lo contengono(2)
	map<unsigned int, vector<unsigned int>> faccePerSpigolo;

	// Itera su tutte le facce
	for (size_t i = 0; i < Originale.SpigoliFacce.size(); ++i) {
		for (auto id_spigolo : Originale.SpigoliFacce[i]) {
			faccePerSpigolo[id_spigolo].push_back(i);
		}
	}

	//nuove strutture dati
	unsigned int numSpigoli = faccePerSpigolo.size();   // Numero di spigoli duali = numero di spigoli originali
	Duale.IdSpigoli.resize(numSpigoli);
	Duale.EstremiSpigoli = MatrixXi(2, numSpigoli);

	// Itera sulla mappa per salvare i nuovi spigoli (estremi = facce adiacenti nel poliedro originale)
	unsigned int idx = 0;
	for (const auto& coppia : faccePerSpigolo) {
		Duale.IdSpigoli[idx] = coppia.first;  // stesso id dello spigolo originale
   
		// Le due facce adiacenti diventano gli estremi dello spigolo nel duale
		if (coppia.second.size() == 2) {
			Duale.EstremiSpigoli(0, idx) = coppia.second[0];
			Duale.EstremiSpigoli(1, idx) = coppia.second[1];
		} else {
			// ERRORE NUMERO FACCE ADIACENTI
			cout << "piu di 2 facce adiacenti" << endl;
		}
		idx++;
	}

	// FACCE

	// Mappa: vertice originale a facce che lo contengono
	map<unsigned int, vector<unsigned int>> faccePerVertice;

	for (size_t i = 0; i < Originale.VerticiFacce.size(); ++i) {
		for (auto id_vertice : Originale.VerticiFacce[i]) {
			faccePerVertice[id_vertice].push_back(i);  // la faccia i contiene il vertice
		}
	}

	// Crea le facce del duale
	unsigned int numFacceDual = faccePerVertice.size();
	Duale.IdFacce.resize(numFacceDual);
	Duale.VerticiFacce.resize(numFacceDual);
	Duale.SpecificaFacce.resize(numFacceDual);

	idx = 0;
	for (const auto& [id_vertice, facce] : faccePerVertice) {
		Duale.IdFacce[idx] = idx;  // nuovo id faccia duale
		Duale.VerticiFacce[idx] = facce;  // ogni faccia è composta dai baricentri (ID facce originali)
		idx++;
	}


    return Duale;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned int aggiunta_vertice(MatrixXd& CoordinateVertici, vector<unsigned int>& IdVertici, unsigned int& contatore_idvertici, const Vector3d& coord)
{
    for (unsigned int i = 0; i < CoordinateVertici.cols(); ++i) {
        if ((CoordinateVertici.col(i) - coord).norm() < 1e-9) {
            return i;
        }
    }

    // Aggiunge nuovo vertice
    CoordinateVertici.conservativeResize(3, contatore_idvertici + 1);
    CoordinateVertici.col(contatore_idvertici) = coord;
    IdVertici.push_back(contatore_idvertici);
    return contatore_idvertici++;
}


void aggiunta_spigolo(MatrixXi& EstremiSpigoli, vector<unsigned int>& IdSpigoli, unsigned int& contatore_idspigoli, map<pair<unsigned int, unsigned int>, unsigned int>& spigoli_id_map, unsigned int v1, unsigned int v2)
{
    if (v1 > v2) swap(v1, v2);
    pair<unsigned int, unsigned int> spigolo = {v1, v2};

    if (spigoli_id_map.find(spigolo) == spigoli_id_map.end()) {
        EstremiSpigoli.conservativeResize(2, contatore_idspigoli + 1);
        EstremiSpigoli.col(contatore_idspigoli) << v1, v2;
        IdSpigoli.push_back(contatore_idspigoli);
        spigoli_id_map[spigolo] = contatore_idspigoli;
        ++contatore_idspigoli;
    }
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
// funzione per triangolare le facce di un poliedro
Poliedro Triangolato(const Poliedro& input, const unsigned int& b) {
    //inizializiazione del poliedro triangolato, dei contatori e di una mappa per gli spigoli  
	Poliedro triangolato;
    unsigned int contatore_idvertici = 0;
    unsigned int contatore_idspigoli = 0;
	unsigned int contatore_idfacce = 0;
    map<pair<unsigned int, unsigned int>, unsigned int> spigoli_id_map;

	// iterazione logica per ogni faccia
    for (const auto& f_idx : input.FaccePoliedro[0]) {
        const auto& face_vertices_indices = input.VerticiFacce[f_idx];
        assert(face_vertices_indices.size() == 3); // controllo che le facce siano triangolari

		// recupero dei vertici per una specifica faccia
        Vector3d A = input.CoordinateVertici.col(face_vertices_indices[0]);
        Vector3d B = input.CoordinateVertici.col(face_vertices_indices[1]);
        Vector3d C = input.CoordinateVertici.col(face_vertices_indices[2]);
		
		// creazione di una griglia per i nuovi vertici
        vector<vector<unsigned int>> griglia_vertici(b + 1, vector<unsigned int>(b + 1));

        // Creazione dei nuovi vertici sulla faccia
        for (int i = 0; i <= b; ++i) {
            for (int j = 0; j <= b - i; ++j) {
                double u = static_cast<double>(i) / b;
                double v = static_cast<double>(j) / b;
                double w = static_cast<double>(b - i - j) / b;
                Vector3d P = u * A + v * B + w * C;
                griglia_vertici[i][j] = aggiunta_vertice(triangolato.CoordinateVertici, triangolato.IdVertici, contatore_idvertici, P);
            }
        }

        // Creazione delle nuove facce triangolari e degli spigoli
        for (int i = 0; i < b; ++i) {
            for (int j = 0; j < b - i; ++j) {
                unsigned int v0 = griglia_vertici[i][j];
                unsigned int v1 = griglia_vertici[i + 1][j];
                unsigned int v2 = griglia_vertici[i][j + 1];

                triangolato.IdFacce.push_back(contatore_idfacce++);
                triangolato.VerticiFacce.push_back({v0, v1, v2});
				
                aggiunta_spigolo(triangolato.EstremiSpigoli, triangolato.IdSpigoli, contatore_idspigoli, spigoli_id_map, v0, v1);
				aggiunta_spigolo(triangolato.EstremiSpigoli, triangolato.IdSpigoli, contatore_idspigoli, spigoli_id_map, v1, v2);
				aggiunta_spigolo(triangolato.EstremiSpigoli, triangolato.IdSpigoli, contatore_idspigoli, spigoli_id_map, v2, v0);
                

                if (i + j + 1 < b) {
                    unsigned int v3 = griglia_vertici[i + 1][j + 1];
                    triangolato.IdFacce.push_back(contatore_idfacce++);
                    triangolato.VerticiFacce.push_back({v1, v3, v2});
                    aggiunta_spigolo(triangolato.EstremiSpigoli, triangolato.IdSpigoli, contatore_idspigoli, spigoli_id_map, v1, v3);
                    aggiunta_spigolo(triangolato.EstremiSpigoli, triangolato.IdSpigoli, contatore_idspigoli, spigoli_id_map, v3, v2);
                    aggiunta_spigolo(triangolato.EstremiSpigoli, triangolato.IdSpigoli, contatore_idspigoli, spigoli_id_map, v2, v1);
                }
            }
        }
    }
	
	// completamento restante parte della struttura dati 
    triangolato.IdPoliedro = {0};
    triangolato.FaccePoliedro = {{0, contatore_idfacce}};
    triangolato.NumVerticiPoliedro = triangolato.IdVertici.size();
    triangolato.NumSpigoliPoliedro = triangolato.IdSpigoli.size();
    triangolato.NumFaccePoliedro = contatore_idfacce;
    triangolato.VerticiPoliedro = {triangolato.IdVertici};
    triangolato.SpigoliPoliedro = {triangolato.IdSpigoli};
	PolopolaSpigoliFacce(triangolato);

    return triangolato;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PolopolaSpigoliFacce(Poliedro& poliedro) {
    // inizializzazione di SpigoliFacce
    poliedro.SpigoliFacce.clear();
    poliedro.SpigoliFacce.resize(poliedro.IdFacce.size()); // Pre-alloca spazio per tutte le facce

    // creazione di una mappa per la ricerca efficiente degli ID degli spigoli
    map<pair<unsigned int, unsigned int>, unsigned int> vertici_to_spigolo_id;

    for (unsigned int i = 0; i < poliedro.IdSpigoli.size(); ++i) {
        unsigned int v1 = poliedro.EstremiSpigoli(0, i);
        unsigned int v2 = poliedro.EstremiSpigoli(1, i);

        // normalizza la coppia: il vertice con ID minore viene prima
        if (v1 > v2) {
            std::swap(v1, v2);
        }
        vertici_to_spigolo_id[{v1, v2}] = poliedro.IdSpigoli[i]; // memorizza l'ID dello spigolo
    }

    // iterazione su ogni faccia per trovare i suoi spigoli
    for (unsigned int f_idx = 0; f_idx < poliedro.IdFacce.size(); ++f_idx) {
        const auto& vertici_faccia = poliedro.VerticiFacce[f_idx]; // ottenimento dei vertici della faccia corrente
        unsigned int num_vertici_faccia = vertici_faccia.size();

        // ogni faccia è un poligono, quindi gli spigoli sono formati da coppie di vertici consecutivi più lo spigolo che collega l'ultimo vertice al primo
        for (unsigned int i = 0; i < num_vertici_faccia; ++i) {
            unsigned int v_curr = vertici_faccia[i];
            unsigned int v_next = vertici_faccia[(i + 1) % num_vertici_faccia]; // per chiudere il ciclo (ultimo con primo)

            if (v_curr > v_next) {
                swap(v_curr, v_next);
            }

            // ricerca dell'ID dello spigolo nella mappa
            auto it = vertici_to_spigolo_id.find({v_curr, v_next});
            if (it != vertici_to_spigolo_id.end()) {
                // Spigolo trovato, aggiungi il suo ID alla lista degli spigoli di questa faccia
                poliedro.SpigoliFacce[f_idx].push_back(it->second);
            } else {
                // errore in caso di un'inconsistenza nel poliedro
                cerr << "Attenzione: Spigolo tra vertici " << v_curr << " e " << v_next
                          << " per la faccia " << f_idx << " non trovato nella mappa degli spigoli.\n";
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// funzione per esportare principali caratteristiche del poliedro
void EsportaCelle(const Poliedro& poliedro) {
	ofstream Celle0Ds("Cells0Ds.txt");
	ofstream Celle1Ds("Cells1Ds.txt");
	ofstream Celle2Ds("Cells2Ds.txt");

    Celle0Ds << "Id:x;y;z\n";
    for (size_t i = 0; i < poliedro.IdVertici.size(); ++i) {
        unsigned int id = poliedro.IdVertici[i];
        Celle0Ds << id << ":"
                 << setprecision(16)
                 << poliedro.CoordinateVertici(0, id) << ";"
                 << poliedro.CoordinateVertici(1, id) << ";"
                 << poliedro.CoordinateVertici(2, id) << "\n";
    }

    Celle1Ds << "Id:inizio;fine\n";
    for (size_t i = 0; i < poliedro.IdSpigoli.size(); ++i) {
        unsigned int id = poliedro.IdSpigoli[i];
        Celle1Ds << id << ":"
                 << poliedro.EstremiSpigoli(0, id) << ";"
                 << poliedro.EstremiSpigoli(1, id) << "\n";
    }

    Celle2Ds << "Id:vertici\n";
    for (size_t i = 0; i < poliedro.IdFacce.size(); ++i) {
        unsigned int id = poliedro.IdFacce[i];
        Celle2Ds << id << ":";
        for (unsigned int v : poliedro.VerticiFacce[i]) {
            Celle2Ds << ";" << v;
        }
        Celle2Ds << "\n";
    }

    Celle0Ds.close();
    Celle1Ds.close();
    Celle2Ds.close();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// funzione per creare un grafo di adiacenza del poliedro
vector<vector<pair<unsigned int, double>>> Adiacenza(const Poliedro& P) {

    const unsigned int numVertici = P.IdVertici.size();
    const unsigned int numSpigoli = P.IdSpigoli.size();

    vector<vector<pair<unsigned int, double>>> adiacenza(numVertici); // alla posizione j ci sono tutti i vertici adiacenti al vertice j con relativa distanza

    for (unsigned int i = 0; i < numSpigoli; ++i) {
        unsigned int a = P.EstremiSpigoli(0, i);
        unsigned int b = P.EstremiSpigoli(1, i);

        Vector3d va = P.CoordinateVertici.col(a);
        Vector3d vb = P.CoordinateVertici.col(b);

        double distanza = (va - vb).norm();

        adiacenza[a].emplace_back(b, distanza);
        adiacenza[b].emplace_back(a, distanza);
    }

    return adiacenza;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// funzione che calcola il percorso minimo fra due vertici del poliedro 
void PercorsoMinimo(Poliedro& P, unsigned int& origine, unsigned int& estremo, vector<double>& bufferSpigoli) {

    const unsigned int N = P.IdVertici.size(); // numero di vertici totali 
    const auto& adiacenza = Adiacenza(P);  // grafo di adiancenza

    vector<double> distanza(N, numeric_limits<double>::infinity()); // dichiarazione vettore distanza con N componenti inizializzate a infinito
    vector<int> predecessore(N, -1); // dichiarazione del vettore dei predecessori (inizializzati a -1; termine di non definizione)
    priority_queue<pair<double, unsigned int>, vector<pair<double, unsigned int>>, greater<>> coda; // dichiarazione della coda di priorità

    distanza[origine] = 0.0; // inizializzazione della distanza a 0
    coda.emplace(0.0, origine); 

    // implementazione algoritmo di Dijkstra per la ricerca del percorso minimo fra due nodi
    while (!coda.empty()) {
        auto [dist_corrente, nodo] = coda.top();
        coda.pop();

        if (nodo == estremo) break;

        for (const auto& [vicino, peso] : adiacenza[nodo]) {
            double nuova_distanza = dist_corrente + peso;
            if (nuova_distanza < distanza[vicino]) {
                distanza[vicino] = nuova_distanza;
                predecessore[vicino] = nodo;
                coda.emplace(nuova_distanza, vicino);
            }
        }
    }

    // Ricostruzione percorso minimo
    vector<unsigned int> percorso;
    for (int v = estremo; v != -1; v = predecessore[v]) {
        percorso.push_back(v);
    }
    reverse(percorso.begin(), percorso.end());

    // inizializzazione buffer come vettori di 0 (colore blu)
    bufferSpigoli.assign(P.IdSpigoli.size(), 0.0);

    // inizializzazione variabili per la misura della lunghezza del percorso
    double lunghezza_totale = 0.0;
    int num_archi = 0;
	
	// calcolo della lunghezza del percorso e assegnazione del valore 2 (colore rosso) agli spigoli che fanno parte del percorso 
    for (size_t i = 1; i < percorso.size(); ++i) {
        unsigned int u = percorso[i - 1];
        unsigned int v = percorso[i];

        for (size_t id = 0; id < P.IdSpigoli.size(); ++id) {
            auto a = P.EstremiSpigoli(0, id);
            auto b = P.EstremiSpigoli(1, id);
            if ((a == u && b == v) || (a == v && b == u)) {
                bufferSpigoli[id] = 2.0;
                // Calcolo distanza euclidea tra i due vertici
                Vector3d p1 = P.CoordinateVertici.col(a);
                Vector3d p2 = P.CoordinateVertici.col(b);
                lunghezza_totale += (p1 - p2).norm();
                num_archi++;
                break;
            }
        }
    }

    // stampa a terminale
    cout << "Numero di archi nel percorso minimo: " << num_archi << endl;
    cout << "Lunghezza totale del percorso minimo: " << lunghezza_totale << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// funzione che compila il vettore di UCDProperty sulla base dei dati salavati nel buffer
void Colorazione(const std::vector<double>& bufferSpigoli, std::vector<UCDProperty<double>>& ColoreSpigoli) {
    ColoreSpigoli.clear();

    ColoreSpigoli.push_back(UCDProperty<double>{
        "Colorazione Percorso Minimo",  // Nome proprietà
        "",
        static_cast<unsigned int>(bufferSpigoli.size()),
        1,                     // Scalare per spigolo
        bufferSpigoli.data()
    });
}

}
