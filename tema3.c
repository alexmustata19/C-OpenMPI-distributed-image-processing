/*
Mustata Alexandru-Ionut
Grupa 331CB
Tema 3 APD
Prelucrare de imagini folosind retele MPI
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <unistd.h>

/*TAG-urile necesare pentru marcarea mesajelor.*/
#define TAG_TOPOLOGIE 0
#define TAG_NUMAR_DE_LINII 1
#define TAG_NUMAR_DE_COLOANE 2
#define TAG_FILTRU 3
#define TAG_SHUTDOWN 4
#define TAG_DATE 5
#define TAG_STATISTICA 6
#define TAG_MAX_PIXEL_VALUE 7

/*Etichetele pentru cele 4 tipuri de filtre si
declararea lor.*/
#define FILTRU_SMOOTH 10
#define FILTRU_BLUR 11
#define FILTRU_SHARPEN 12
#define FILTRU_MEAN_REMOVAL 13

int smooth[] = {
	1, 1, 1,
	1, 1, 1,
	1, 1, 1
};
int smoothS = 9;

int blur[] = {
	1, 2, 1,
	2, 4, 2,
	1, 2, 1
};
int blurS = 16; 

int sharpen[] = {
	0, -2, 0,
	-2, 11, -2,
	0, -2, 0
};
int sharpenS = 3;

int mean_removal[] = {
	-1, -1, -1,
	-1, 9, -1,
	-1, -1, -1
};
int mean_removalS = 1;

/*Structura ce reprezinta formatul unei imagini PGM.*/
typedef struct image{
	char* antet;
	char* linieComentariu;
	int nrLinii, nrColoane, maxPixelValue;
	int** matrice;
}Timage;

/*Functia realizeaza citirea din fisier a unei imagini PGM.*/
Timage* citireImagine(char* numeImagineIn){
	/*Citirea liniei cu eticheta P2.*/
	FILE* f = fopen(numeImagineIn, "r");
	size_t l;
	char* linie = (char*)malloc(10000 * sizeof(int));
	getline(&linie, &l, f);
	if(strcmp(linie,"P2\n") != 0){
		return NULL;
	}

	/*Citirea eventualei linii de comentariu si a numarului
	de linii si coloane si a valorii maxime a pixelilor.*/
	Timage* img = (Timage*)malloc(sizeof(Timage));
	img->antet = (char*)malloc(strlen(linie)*sizeof(char));
	strcpy(img->antet, linie);

	getline(&linie, &l, f);
	if(linie[0] == '#'){
		img->linieComentariu = (char*)malloc(strlen(linie)*sizeof(char));
		strcpy(img->linieComentariu, linie);
		fscanf(f, "%d", &(img->nrColoane));
	}else{
		img->linieComentariu = (char*)calloc(strlen("nimic") + 1, sizeof(char));
		strcpy(img->linieComentariu, "nimic");
		sscanf(linie, "%d", &(img->nrColoane));
	}
	fscanf(f, "%d", &(img->nrLinii));
	fscanf(f, "%d", &(img->maxPixelValue));
	img->maxPixelValue = 255;

	/*Citirea matricii ce reprezinta imaginea.*/
	img->matrice = (int**)malloc((img->nrLinii+2) * sizeof(int*));
	int i = 0;
	for(i = 0; i < img->nrLinii + 2; i++){
		img->matrice[i] = (int*)malloc((img->nrColoane + 2) * sizeof(int));
	}
	int j = 0;
	for(i = 1; i < img->nrLinii + 1; i++){
		for(j = 1; j < img->nrColoane + 1; j++){
			fscanf(f, "%d", &(img->matrice[i][j]));
		}
	}
	for(i = 0; i < img->nrLinii + 2; i++){
		img->matrice[i][0] = img->matrice[i][img->nrColoane+1] = 0;
	}
	for(i = 0; i < img->nrColoane + 2; i++){
		img->matrice[0][i] = img->matrice[img->nrLinii+1][i] = 0;
	}

	fclose(f);
	free(linie);
	return img;
}

/*Functia scrie o imagine PGM in fisier.*/
void scriereImagine(char* numeImagineOut, Timage* img){
	FILE* f = fopen(numeImagineOut, "w");
	fprintf(f, "%s", img->antet);
	if(img->linieComentariu[0] == '#'){
		fprintf(f, "%s", img->linieComentariu);
	}
	fprintf(f, "%d ", img->nrColoane);
	fprintf(f, "%d\n", img->nrLinii);
	fprintf(f, "%d\n", img->maxPixelValue);
	int i = 0;
	int j = 0;
	for(i = 1; i < img->nrLinii+1; i++){
		for(j = 1; j < img->nrColoane+1; j++){
			fprintf(f, "%d\n", img->matrice[i][j]);
		}
	}
	fclose(f);
}

/*Functia dezaloca memoria alocata pentru o imagine.*/
void freeImagine(Timage** img){
	Timage* aux = *img;
	free(aux->antet);
	free(aux->linieComentariu);
	int i = 0;
	for(i = 0; i < aux->nrLinii+2; i++){
		free(aux->matrice[i]);
	}
	free(aux->matrice);
	free(aux);
	*img = NULL;
}

/*Functia realizeaza stabilirea topologiei care va fi folosita la modificarea imaginilor.
In vectorul copii, pozitiile ce vor fi ocupate de 1 vor reprezenta faptul ca acel nod
este copilul procesului ce a apelat functia. Se returneaza parintele nodului ce a
apelat functia.*/
int stabilireTopologie(char* numeFisier, int rank, int* copii, int nProcese){
	int rootConected = 0;
	int parinte = -1;
	if(rank == 0){
		rootConected = 1;
	}

	int* vecini = (int*)malloc(nProcese * sizeof(int));
	int i = 0;
	for(i = 0; i < nProcese; i++){
		vecini[i] = 0;
	}

	/*Citire topologie.*/
	FILE* f = fopen(numeFisier, "r");
	char* linie = (char*)malloc(10000 * sizeof(int));
	int buffer = -1;

	/*Gasirea liniei corespunzatoare procesului.*/
	while(buffer != rank){
		size_t l;
		getline(&linie, &l, f);
		sscanf(linie, "%d", &buffer);
	}

	fclose(f);

	/*Citirea vecinilor de pe linie.*/
	char* rest;
	char* numar;
	strtok_r(linie, " ", &rest);
	while(numar = strtok_r(rest, " ", &rest)){
		buffer = atoi(numar);
		vecini[buffer] = 1;
	}

	/*Comunicarea cu vecinii in vederea stabilirii topologiei.*/
	buffer = -1;
	for(i = 0; i < nProcese; i++){
		int j = 0;
		for(j = 0; j < nProcese; j++){
			if(vecini[j] == 1){
				MPI_Sendrecv(&rootConected, 1, MPI_INT, j, TAG_TOPOLOGIE, &buffer, 1, MPI_INT, j, TAG_TOPOLOGIE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if(rootConected == 0 && buffer == 1){
					/*Daca vecinul ajunge la radacina si procesul curent nu,
					atunci se seteaza vecinul ca parinte, iar procesul curent
					va ajunge prin el la radacina.*/
					rootConected = 1;
					parinte = j;
				}else if(rootConected == 1 && buffer == 0){
					/*Daca vecinul nu ajunge la radacina si procesul curent da,
					atunci vecinul devine copilul procesului curent.*/
					copii[j] = 1;
				}
			}
		}
	}
	free(vecini);
	free(linie);
	return parinte;
}

/*Functia aplica filtrul specificat pe o matrice de nrLinii X nrColoane obtinand 
calculul in matricea rezultat. Elementele matrici vor avea in urma calcului valori
in [0, maxPixelValue].*/
void aplicareFiltru(int** matrice, int** rezultat, int nrLinii, int nrColoane, int tipFiltru, int maxPixelValue){
	int* filtru = NULL;
	int filtruS = 0;
	if(tipFiltru == FILTRU_SMOOTH){
		filtru = smooth;
		filtruS = smoothS;
	}else if(tipFiltru == FILTRU_BLUR){
		filtru = blur;
		filtruS = blurS;
	}else if(tipFiltru == FILTRU_SHARPEN){
		filtru = sharpen;
		filtruS = sharpenS;
	}else if(tipFiltru == FILTRU_MEAN_REMOVAL){
		filtru = mean_removal;
		filtruS = mean_removalS;
	}
	int i = 0;
	int j = 0;
	for(i = 1; i < nrLinii + 1; i++){
		for(j = 1; j < nrColoane + 1; j++){
			float aux = 1.0 * (filtru[0] * matrice[i-1][j-1] + filtru[1] * matrice[i-1][j] + filtru[2] * matrice[i-1][j+1] +
				filtru[3] * matrice[i][j-1] + filtru[4] * matrice[i][j] + filtru[5] * matrice[i][j+1] +
				filtru[6] * matrice[i+1][j-1] + filtru[7] * matrice[i+1][j] + filtru[8] * matrice[i+1][j+1]) / filtruS;
			if(aux <= 0){
				rezultat[i][j] = 0;
			}else if(aux >= maxPixelValue){
				rezultat[i][j] = maxPixelValue;
			}else{
				rezultat[i][j] = (int)aux;
			}
		}
	}
}

int main(int argc, char **argv){
	int rank;
	int nProcese;

	if(argc != 4){
		printf("Numar incorect de argumente!\n");
		return 0;
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcese);

	/*Stabilirea topologiei*/
	int* copii = (int*)malloc(nProcese * sizeof(int));
	int i = 0;
	for(i = 0; i < nProcese; i++){
		copii[i] = 0;
	}
	int parinte = stabilireTopologie(argv[1], rank, copii, nProcese);

	/*Calcularea numarului de copii si memorarea indicilor
	acestora intr-un vector.*/
	int nrCopii = 0;
	for(i = 0; i < nProcese; i++){
		nrCopii += copii[i];
	}
	int* indiciCopii = NULL;
	if(nrCopii > 0){
		indiciCopii = (int*)malloc(nrCopii * sizeof(int));
		int index = 0;
		for(i = 0; i < nProcese; i++){
			if(copii[i] == 1){
				indiciCopii[index] = i;
				index++;
			}
		}
	}

	/*Executia programului in functie de tipul de nod din retea.*/
	if(rank == 0){
		/*Logica radacinii. Citirea imaginilor din fisierul cu lista
		si trimiterea acestora spre procesare in retea. La final se
		trimite semnalul de SHUTDOWN si se colecteaza statistica.*/
		
		/*Determinare numarului de task-uri.*/
		FILE* f = fopen(argv[2],"r");
		int taskuri  = 0;
		char c;
		fscanf(f, "%d%c", &taskuri, &c);

		for(i = 0; i < taskuri; i++){
			/*Citirea imaginii aferente fiecarui task si
			stabilirea tipului de filtru.*/
			size_t l;
			char* linie = (char*)malloc(10000 * sizeof(int));
			getline(&linie, &l, f);
			char* rest;
			char* filtru = strtok_r(linie, " \n\t", &rest);
			int tipFiltru = 0;
			if(strcmp(filtru, "smooth") == 0){
				tipFiltru = FILTRU_SMOOTH;
			}else if(strcmp(filtru, "blur") == 0){
				tipFiltru = FILTRU_BLUR;
			}else if(strcmp(filtru, "sharpen") == 0){
				tipFiltru = FILTRU_SHARPEN;
			}else if(strcmp(filtru, "mean_removal") == 0){
				tipFiltru = FILTRU_MEAN_REMOVAL;
			}
			char* numeImagineIn = strtok_r(rest, " \n\t", &rest);
			char* numeImagineOut = strtok_r(rest, " \n\t", &rest);
			Timage* img = citireImagine(numeImagineIn);
			if(img == NULL){
				exit(0);
			}
			int nrLinii = img->nrLinii;

			/*Impartirea imagini de procesat la nodurile copii.*/
			if(nrLinii >= nrCopii){
				/*Daca exista mai multe linii decat copii atunci primii
				nrCopii-1 copii vor primi [nrLinii/nrCopii], iar ultimul
				copil va primi restul de linii.*/
				float aux = 1.0 * img->nrLinii / nrCopii;
				int felie = floor(aux);
				int j = 0;
				int k = 0;
				int begin = 1;
				int end = felie+1;

				/*Trimiterea bucatilor egale la primii copii.*/
				for(j = 0; j < nrCopii-1; j++){
					MPI_Send(&felie, 1, MPI_INT, indiciCopii[j], TAG_NUMAR_DE_LINII, MPI_COMM_WORLD);
					MPI_Send(&(img->nrColoane), 1, MPI_INT, indiciCopii[j], TAG_NUMAR_DE_COLOANE, MPI_COMM_WORLD);
					MPI_Send(&tipFiltru, 1, MPI_INT, indiciCopii[j], TAG_FILTRU, MPI_COMM_WORLD);
					MPI_Send(&(img->maxPixelValue), 1, MPI_INT, indiciCopii[j], TAG_MAX_PIXEL_VALUE, MPI_COMM_WORLD);
					for(k = begin - 1; k < end + 1; k++){
						MPI_Send(img->matrice[k], img->nrColoane+2, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD);
					}
					begin += felie;
					end += felie;
				}
				
				/*Trimiterea ultimei bucati rest la ultimul copil.*/
				int felieRest = img->nrLinii - (felie * (nrCopii-1));
				MPI_Send(&felieRest, 1, MPI_INT, indiciCopii[nrCopii-1], TAG_NUMAR_DE_LINII, MPI_COMM_WORLD);
				MPI_Send(&(img->nrColoane), 1, MPI_INT, indiciCopii[nrCopii-1], TAG_NUMAR_DE_COLOANE, MPI_COMM_WORLD);
				MPI_Send(&tipFiltru, 1, MPI_INT, indiciCopii[nrCopii-1], TAG_FILTRU, MPI_COMM_WORLD);
				MPI_Send(&(img->maxPixelValue), 1, MPI_INT, indiciCopii[nrCopii-1], TAG_MAX_PIXEL_VALUE, MPI_COMM_WORLD);
				begin = 1 + (nrCopii-1) * felie;
				end = nrLinii + 1;
				for(k = begin - 1; k < end + 1; k++){
					MPI_Send(img->matrice[k], img->nrColoane+2, MPI_INT, indiciCopii[nrCopii-1], TAG_DATE, MPI_COMM_WORLD);
				}

				/*Primirea rezultatelor de la primii copii.*/
				int count = 1;
				for(j = 0; j < nrCopii-1; j++){
					for(k = 0; k < felie; k++){
						MPI_Recv(img->matrice[count] + 1, img->nrColoane, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						count++;
					}
				}

				/*Primirea rezultatului de la ultimul copil.*/
				for(k = 0; k < felieRest; k++){
					MPI_Recv(img->matrice[count] + 1, img->nrColoane, MPI_INT, indiciCopii[nrCopii-1], TAG_DATE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					count++;
				}
			}else{
				/*Numarul de linii este mai mic decat numarul de copii,
				deci primii nrLinii copii primesc cate o linie, restul
				vor primi zero pentru a sti ca la acest task nu lucreaza.*/
				int j = 0;

				/*Trimiterea cate unei linii la primii copii.*/
				int unitate = 1;
				for(j = 0; j < nrLinii; j++){
					MPI_Send(&unitate, 1, MPI_INT, indiciCopii[j], TAG_NUMAR_DE_LINII, MPI_COMM_WORLD);
					MPI_Send(&(img->nrColoane), 1, MPI_INT, indiciCopii[j], TAG_NUMAR_DE_COLOANE, MPI_COMM_WORLD);
					MPI_Send(&tipFiltru, 1, MPI_INT, indiciCopii[j], TAG_FILTRU, MPI_COMM_WORLD);
					MPI_Send(&(img->maxPixelValue), 1, MPI_INT, indiciCopii[j], TAG_MAX_PIXEL_VALUE, MPI_COMM_WORLD);
					MPI_Send(img->matrice[j], img->nrColoane+2, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD);
					MPI_Send(img->matrice[j+1], img->nrColoane+2, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD);
					MPI_Send(img->matrice[j+2], img->nrColoane+2, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD);
				}
				
				/*Trimiterea instiintarii de neparticipare la task.*/
				unitate = 0;
				for(j = nrLinii; j < nrCopii; j++){
					MPI_Send(&unitate, 1, MPI_INT, indiciCopii[j], TAG_NUMAR_DE_LINII, MPI_COMM_WORLD);
				}

				/*Primirea liniilor procesate de la primii copii.*/
				for(j = 0; j < nrLinii; j++){
					MPI_Recv(img->matrice[j+1], img->nrColoane+2, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			}

			/*Scrierea imaginii procesate in fisierul indicat de output.*/
			scriereImagine(numeImagineOut, img);
			freeImagine(&img);
			free(linie);
		}

		/*Dupa terminarea task-urilor, se declanseaza etapa de SHUTDOWN.
		Se trimite in retea semnalul, iar ca raspuns de la fiecare copil
		se asteapta statistica din subarborele dominat de el sub forma
		unui vector.*/

		/*Initierea semnalului de SHUTDOWN.*/
		for(i = 0; i < nrCopii; i++){
			MPI_Send(&rank, 1, MPI_INT, indiciCopii[i], TAG_SHUTDOWN, MPI_COMM_WORLD);
		}

		/*Primirea si agregarea statisticilor de la copii
		intr-un singur vector.*/
		int* vec = (int*)malloc(nProcese * sizeof(int));
		for(i = 0; i < nProcese; i++){
			vec[i] = -1;
		}
		vec[rank] = 0;
		int j = 0;
		int* vecBuf = (int*)malloc(nProcese * sizeof(int));
		for(i = 0; i < nrCopii; i++){
			MPI_Recv(vecBuf, nProcese, MPI_INT, indiciCopii[i], TAG_STATISTICA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(j = 0; j < nProcese; j++){
				if(vecBuf[j] != -1){
					vec[j] = vecBuf[j];
				}
			}
		}

		/*Scrierea statistici in fisierul corespunzator.*/
		FILE* stats = fopen(argv[3], "w");
		for(i = 0; i < nProcese; i++){
			fprintf(stats, "%d: %d\n", i, vec[i]);
		}
		free(vec);
		free(vecBuf);
		fclose(stats);
		fclose(f);
	}else if(nrCopii != 0){
		/*Logica nodulurilor intermediare. Un nod de acest tip primeste de la
		parinte o bucata de imagine care o imparte la copii pentru procesare.
		Rezultatele individuale primite de la copii sunt grupate si trimise
		parintelui ca raspuns.*/
		while(1){
			int nrLinii;
			int nrColoane;
			int tipFiltru;
			int maxPixelValue;
			MPI_Status status;
			MPI_Recv(&nrLinii, 1, MPI_INT, parinte, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if(status.MPI_TAG == TAG_SHUTDOWN){
				/*Daca mesajul primit are tag-ul de SHUTDOWN, se transmite copiilor acest mesaj
				si se asteapta ca raspuns de la fiecare un vector pentru statistica. Dupa aceea,
				se trimite vectorul rezultat de statistica parintelui si se termina executia
				procesului.*/
				for(i = 0; i < nrCopii; i++){
					MPI_Send(&nrLinii, 1, MPI_INT, indiciCopii[i], TAG_SHUTDOWN, MPI_COMM_WORLD);
				}
				int* vec = (int*)malloc(nProcese * sizeof(int));
				for(i = 0; i < nProcese; i++){
					vec[i] = -1;
				}
				vec[rank] = 0;
				int j = 0;
				int* vecBuf = (int*)malloc(nProcese * sizeof(int));
				for(i = 0; i < nrCopii; i++){
					MPI_Recv(vecBuf, nProcese, MPI_INT, indiciCopii[i], TAG_STATISTICA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for(j = 0; j < nProcese; j++){
						if(vecBuf[j] != -1){
							vec[j] = vecBuf[j];
						}
					}
				}
				MPI_Send(vec, nProcese, MPI_INT, parinte, TAG_STATISTICA, MPI_COMM_WORLD);
				free(vec);
				free(vecBuf);
				break;
			}else if(nrLinii == 0){
				/*Daca nu se primesc linii, se va trece la asteptarea urmatorului
				mesaj care va sosi de la parinte.*/
				continue;
			}else{
				/*Trebuie receptionat si numarul de coloane, tipul filtrului, valoarea maxima a pixelilor
				din imaginine si bucata efectiva din imagine ce va fi impartita copiilor.*/
				MPI_Recv(&nrColoane, 1, MPI_INT, parinte, TAG_NUMAR_DE_COLOANE, MPI_COMM_WORLD, &status);
				MPI_Recv(&tipFiltru, 1, MPI_INT, parinte, TAG_FILTRU, MPI_COMM_WORLD, &status);
				MPI_Recv(&maxPixelValue, 1, MPI_INT, parinte, TAG_MAX_PIXEL_VALUE, MPI_COMM_WORLD, &status);
				int** matrice = (int**)malloc((nrLinii+2) * sizeof(int*));
				for(i = 0; i < nrLinii+2; i++){
					matrice[i] = (int*)malloc((nrColoane+2) * sizeof(int));
					MPI_Recv(matrice[i], nrColoane+2, MPI_INT, parinte, TAG_DATE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}

				if(nrLinii >= nrCopii){
					/*Daca exista mai multe linii decat copii atunci primii
					nrCopii-1 copii vor primi [nrLinii/nrCopii], iar ultimul
					copil va primi restul de linii.*/
					float aux = 1.0 * nrLinii / nrCopii;
					int felie = floor(aux);
					int j = 0;
					int k = 0;
					int begin = 1;
					int end = felie+1;

					/*Trimiterea bucatilor egale la primii copii.*/
					for(j = 0; j < nrCopii-1; j++){
						MPI_Send(&felie, 1, MPI_INT, indiciCopii[j], TAG_NUMAR_DE_LINII, MPI_COMM_WORLD);
						MPI_Send(&nrColoane, 1, MPI_INT, indiciCopii[j], TAG_NUMAR_DE_COLOANE, MPI_COMM_WORLD);
						MPI_Send(&tipFiltru, 1, MPI_INT, indiciCopii[j], TAG_FILTRU, MPI_COMM_WORLD);
						MPI_Send(&maxPixelValue, 1, MPI_INT, indiciCopii[j], TAG_MAX_PIXEL_VALUE, MPI_COMM_WORLD);
						for(k = begin - 1; k < end + 1; k++){
							MPI_Send(matrice[k], nrColoane+2, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD);
						}
						begin += felie;
						end += felie;
					}

					/*Trimiterea ultimei bucati rest la ultimul copil.*/
					int felieRest = nrLinii - (felie * (nrCopii-1));
					MPI_Send(&felieRest, 1, MPI_INT, indiciCopii[nrCopii-1], TAG_NUMAR_DE_LINII, MPI_COMM_WORLD);
					MPI_Send(&nrColoane, 1, MPI_INT, indiciCopii[nrCopii-1], TAG_NUMAR_DE_COLOANE, MPI_COMM_WORLD);
					MPI_Send(&tipFiltru, 1, MPI_INT, indiciCopii[nrCopii-1], TAG_FILTRU, MPI_COMM_WORLD);
					MPI_Send(&maxPixelValue, 1, MPI_INT, indiciCopii[nrCopii-1], TAG_MAX_PIXEL_VALUE, MPI_COMM_WORLD);
					begin = 1 + (nrCopii-1) * felie;
					end = nrLinii + 1;
					for(k = begin - 1; k < end + 1; k++){
						MPI_Send(matrice[k], nrColoane+2, MPI_INT, indiciCopii[nrCopii-1], TAG_DATE, MPI_COMM_WORLD);
					}

					/*Primirea rezultatelor de la primii copii.*/
					int count = 1;
					for(j = 0; j < nrCopii-1; j++){
						for(k = 0; k < felie; k++){
							MPI_Recv(matrice[count] + 1, nrColoane, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							count++;
						}
					}
					/*Primirea rezultatului de la ultimul copil.*/
					for(k = 0; k < felieRest; k++){
						MPI_Recv(matrice[count] + 1, nrColoane, MPI_INT, indiciCopii[nrCopii-1], TAG_DATE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						count++;
					}
				}else{
					/*Numarul de linii este mai mic decat numarul de copii,
					deci primii nrLinii copii primesc cate o linie, restul
					vor primi zero pentru a sti ca la acest task nu lucreaza.*/
					int j = 0;

					/*Trimiterea cate unei linii la primii copii.*/
					int unitate = 1;
					for(j = 0; j < nrLinii; j++){
						MPI_Send(&unitate, 1, MPI_INT, indiciCopii[j], TAG_NUMAR_DE_LINII, MPI_COMM_WORLD);
						MPI_Send(&nrColoane, 1, MPI_INT, indiciCopii[j], TAG_NUMAR_DE_COLOANE, MPI_COMM_WORLD);
						MPI_Send(&tipFiltru, 1, MPI_INT, indiciCopii[j], TAG_FILTRU, MPI_COMM_WORLD);
						MPI_Send(&maxPixelValue, 1, MPI_INT, indiciCopii[j], TAG_MAX_PIXEL_VALUE, MPI_COMM_WORLD);
						MPI_Send(matrice[j], nrColoane+2, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD);
						MPI_Send(matrice[j+1], nrColoane+2, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD);
						MPI_Send(matrice[j+2], nrColoane+2, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD);
					}

					/*Trimiterea instiintarii de neparticipare la task.*/
					unitate = 0;
					for(j = nrLinii; j < nrCopii; j++){
						MPI_Send(&unitate, 1, MPI_INT, indiciCopii[j], TAG_NUMAR_DE_LINII, MPI_COMM_WORLD);
					}

					/*Primirea liniilor procesate de la primii copii.*/
					for(j = 0; j < nrLinii; j++){
						MPI_Recv(matrice[j+1], nrColoane+2, MPI_INT, indiciCopii[j], TAG_DATE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					}
				}
				/*Trimiterea rezultatului obtinut la parinte.*/
				for(i = 1; i < nrLinii + 1; i++){
					MPI_Send(matrice[i] + 1, nrColoane, MPI_INT, parinte, TAG_DATE, MPI_COMM_WORLD);
				}
				for(i = 0; i < nrLinii+2; i++){
					free(matrice[i]);
				}
				free(matrice);
			}
		}
	}else{
		/*Logica nodulurilor de tip frunza. Un nod de acest tip asteapta
		de la parinte numarul de linii, coloane, tipul filtrului si bucata
		din imagine pe care o va procesa. Daca numarul de linii pe care le
		va primi este zero, inseamna ca la acest task frunza nu va executa
		nicio operatie.*/
		int numarLiniiProcesate = 0;
		while(1){
			int nrLinii;
			int nrColoane;
			int tipFiltru;
			int maxPixelValue;
			MPI_Status status;
			MPI_Recv(&nrLinii, 1, MPI_INT, parinte, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if(status.MPI_TAG == TAG_SHUTDOWN){
				/*Daca se primeste mesajul cu eticheta de SHUTDOWN, atunci frunza
				va crea un vector initializat pe toate pozitiile cu -1. Pe pozitia
				rank-ului ei va pune numarul de linii procesate. Se trimite vectorul
				ca raspuns la semnalul de SHUTDOWN si se termina executia procesului.*/
				int* vec = (int*)malloc(nProcese * sizeof(int));
				for(i = 0; i < nProcese; i++){
					vec[i] = -1;
				}
				vec[rank] = numarLiniiProcesate; 
				MPI_Send(vec, nProcese, MPI_INT, parinte, TAG_STATISTICA, MPI_COMM_WORLD);
				free(vec);
				break;
			}else if(nrLinii == 0){
				/*Daca frunza nu va primi linii de la parinte, atunci va trece la
				asteptarea urmatorului mesaj de la parinte.*/
				continue;
			}else{
				/*Se receptioneaza de la parinte si numarul de coloane, tipul filtrului, valoarea
				maxima a pixelilor din imagine si bucata de imagine de prelucrat.*/
				MPI_Recv(&nrColoane, 1, MPI_INT, parinte, TAG_NUMAR_DE_COLOANE, MPI_COMM_WORLD, &status);
				MPI_Recv(&tipFiltru, 1, MPI_INT, parinte, TAG_FILTRU, MPI_COMM_WORLD, &status);
				MPI_Recv(&maxPixelValue, 1, MPI_INT, parinte, TAG_MAX_PIXEL_VALUE, MPI_COMM_WORLD, &status);
				int** matrice = (int**)malloc((nrLinii+2) * sizeof(int*));
				int** rezultat = (int**)malloc((nrLinii+2) * sizeof(int*));
				for(i = 0; i < nrLinii+2; i++){
					matrice[i] = (int*)malloc((nrColoane+2) * sizeof(int));
					rezultat[i] = (int*)malloc((nrColoane+2) * sizeof(int));
					MPI_Recv(matrice[i], nrColoane+2, MPI_INT, parinte, TAG_DATE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}

				/*Se prelucreaza imaginea si se actualizeaza numarul de linii prelucrate.*/
				aplicareFiltru(matrice, rezultat, nrLinii, nrColoane, tipFiltru, maxPixelValue);
				numarLiniiProcesate += nrLinii;

				/*Se trimite rezultatul obtinut la parinte.*/
				for(i = 1; i < nrLinii + 1; i++){
					MPI_Send(rezultat[i] + 1, nrColoane, MPI_INT, parinte, TAG_DATE, MPI_COMM_WORLD);
				}
				for(i = 0; i < nrLinii+2; i++){
					free(matrice[i]);
					free(rezultat[i]);
				}
				free(matrice);
				free(rezultat);
			}
		}
	}

	/*Dezalocarea vectorilor care au mentinut legaturile
	catre copii.*/
	free(copii);
	if(nrCopii > 0){
		free(indiciCopii);
	}

	MPI_Finalize();
	return 0;
}