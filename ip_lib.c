/*
 * Created by Sebastiano Vascon on 23/03/20.
 * 
 * Gruppo Casteo
 * ID: 69 
 * Simone Jovon  882028
 * Sebastiano Quintavalle  878500
 * Andrea Rosa  882014
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ip_lib.h"
#include "bmp.h"

/**** METODI DI STAMPA ****/

/* Visualizza i dati stampando in ordine le matrici rispetto la terza dimensione */
void ip_mat_show(ip_mat * t) {
	
	if (t) {
		
		unsigned int i, l, j;
		
		printf("Matrix of size %d x %d x %d (hxwxk)\n", t -> h, t -> w, t -> k);
		
		for (l = 0; l < t -> k; l++) {
			printf("Slice %d\n", l);
			for (i = 0; i < t -> h; i++) {
				for (j = 0; j < t -> w; j++) {
					printf("%f ", get_val(t, i, j, l));
				}
				printf("\n");
			}
			printf("\n");
		}
	}
	else {
		printf("Errore ip_mat_show!!!\n");
		exit(1);
	}
}

/* Visualizza a video le statistiche per ogni canale.
 * metodo già implementato */
void ip_mat_show_stats(ip_mat * t){
	
	if (t) {
		
		unsigned int k;
		
		compute_stats(t);
		
		for (k = 0; k < t -> k; k++){
			printf("Channel %d:\n", k);
			printf("\t Min: %f\n",  t -> stat[k].min);
			printf("\t Max: %f\n",  t -> stat[k].max);
			printf("\t Mean: %f\n", t -> stat[k].mean);
		}
	}
	else {
		printf("Errore ip_mat_show_stats!!!\n");
		exit(1);
	}
}

/**** CONVERSIONE BITMAP - IP_MAT ****/

/* Converte una Bitmap in una ip_mat 
 * metodo già implementato */
ip_mat * bitmap_to_ip_mat(Bitmap * img){
	
	if (img) {
		
		unsigned int i=0,j=0;
		
		unsigned char R,G,B;
		
		unsigned int h = img -> h;
		unsigned int w = img -> w;
		
		ip_mat * out = ip_mat_create(h, w, 3, 0);
		
		for (i = 0; i < h; i++)              /* rows */
		{
			for (j = 0; j < w; j++)          /* columns */
			{
				bm_get_pixel(img, j,i,&R, &G, &B);
				set_val(out, i, j, 0, (float) R);
				set_val(out, i, j, 1, (float) G);
				set_val(out, i, j, 2, (float) B);
			}
		}
		
		compute_stats(out);	
		return out;
	}
	else {
		printf("Errore bitmat_to_ip_mat!!!\n");
		exit(1);
	}
}


/* Converte una ip_mat in una bitmap
 * metodo già implementato */
Bitmap * ip_mat_to_bitmap(ip_mat * t) {
	
	if (t) {
		
		Bitmap *b = bm_create(t -> w, t -> h);
		
		unsigned int i, j;
		
		for (i = 0; i < t -> h; i++) {           /* rows */
			for (j = 0; j < t -> w; j++) {       /* columns */
				bm_set_pixel(b, j, i, (unsigned char) get_val(t, i, j, 0),
							 (unsigned char) get_val(t, i, j, 1),
							 (unsigned char) get_val(t, i, j, 2));
			}
		}
		return b;
	}
	else {
		printf("Errore ip_mat_to_bitmap!!!");
		exit(1);
	}
}

/**** RANDOM ****/
/* Genera dei numeri casuali con distribuzione Normale (versione base) */
float get_normal_random(float media, float std){
	
	float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
	float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
	float num = cos(2 * PI * y2) * sqrt(-2. * log(y1));
	
	return media + num * std;
}

/**** PARTE 1: TIPO DI DATI ip_mat E MEMORIA ****/

/* Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali.
 * */
ip_mat * ip_mat_create(unsigned int h, unsigned int w, unsigned  int k, float v) {
	
	/* Contatori */
	unsigned int i, j, c;
	
	/* Allocazione della ip_mat */
	ip_mat *new;
	new = (ip_mat*) malloc(sizeof(ip_mat));
	if (new == NULL) {
		printf("Errore ip_mat_create!!!\n");
		exit(1);
	}
	
	/* Setting delle dimensioni */
	new -> h = h;
	new -> w = w;
	new -> k = k;
	
	/* Allocazione della matrice tridimensionale */
	new -> data = (float***) malloc(sizeof(float**) * h); /* righe */
	if (new -> data == NULL) {
		printf("Errore ip_mat_create!!!\n");
		exit(1);
	}
	for (i = 0; i < h; i++) {
		new -> data[i] = (float**) malloc(sizeof(float*) * w); /* colonne */
		if (new -> data[i] == NULL) {
			printf("Errore ip_mat_create!!!\n");
			exit(1);
		}
		for (j = 0; j < w; j++) {
			new -> data[i][j] = (float*) malloc(sizeof(float) * k); /* canali */
			if (new -> data[i][j] == NULL) {
				printf("Errore ip_mat_create!!!\n");
				exit(1);
			}
			for (c = 0; c < k; c++) {
				set_val(new, i, j, c, v);
			}
		}
	}
	
	/* Allocazione del vettore stat */
	new -> stat = (stats*) malloc(sizeof(stats) * k);
	if (new -> stat == NULL) {
		printf("Errore ip_mat_create!!!\n");
		exit(1);
	}   
	for (c = 0; c < k; c++) {
		new -> stat[c].min  = v;
		new -> stat[c].max  = v;
		new -> stat[c].mean = v;
	}
	
	return new;
}


/* Libera la memoria (data, stat e la struttura) */
void ip_mat_free(ip_mat *a) {
	if (a ){
		
		unsigned int i, j;
		
		/* Deallocazione matrice tridimensionale */
		for (i = 0; i < a -> h; i++) {
			for (j = 0; j < a -> w; j++)
				free(a -> data[i][j]);
			free(a -> data[i]);
		}
		free(a -> data);
		
		/* Deallocazione del vettore stat */
		free(a -> stat);
		
		/* Deallocazione di ip_mat */
		free(a); 
	}
}

/* Restituisce il valore in posizione i, j, k
 * metodo già implementato */
float get_val(ip_mat * a, unsigned int i, unsigned int j, unsigned int k){
	if(a && i < a -> h && j < a -> w && k < a -> k ){  /* j>=0 and k>=0 and i>=0 is non sense*/
		return a -> data[i][j][k];
	}
	else{
		printf("Errore get_val!!!\n");
		exit(1);
	}
}

/* Setta il valore in posizione i, j, k a v
 * metodo già implementato */
void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
	if(a && i < a -> h && j < a -> w && k < a -> k){
		a -> data[i][j][k]=v;
	}
	else{
		printf("Errore set_val!!!\n");
		exit(1);
	}
}

/* Calcola il valore minimo, il massimo e la media per ogni canale
 * e li salva dentro la struttura ip_mat stats */
void compute_stats(ip_mat * t) {
	
	if(t){
		
		unsigned int c;
		
		for (c = 0; c < t  ->  k; c++) { /* scorrimento canali */
			
			unsigned int i, j;
			
			/* Minimo, massimo, accumulatore*/
			float min, max, acc;
			
			/* Inizializzazione delle variabili */
			min = max = get_val (t, 0, 0, c);
			acc = 0.;
			
			for (i = 0; i < t  ->  h; i++) { /* scorrimento righe */
				for (j = 0; j < t  ->  w; j++) { /* scorrimento colonne */
					float temp;
					temp = get_val(t, i, j, c);
					acc += temp;
					if (temp > max)
						max = temp;
					if (temp < min)
						min = temp;
				}
			}
			
			acc /= t  ->  h * t  ->  w; /* media */
			
			/* Setting delle statistiche del canale */
			t  ->  stat[c].min = min;
			t  ->  stat[c].max = max;
			t  ->  stat[c].mean = acc;
			
		}
	}
	else {
		printf("Errore compute_stats!!!\n");
		exit(1);
	}
	
	
}

/* Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento è generato da una gaussiana con media mean e varianza var */
void ip_mat_init_random(ip_mat * t, float mean, float std) {
	
	if(t){
		
		unsigned int i, j, c;
		
		for (c = 0; c < t -> k; c++) { /* scorrimento canali */
			for (i = 0; i < t -> h; i++) { /* scorrimento righe */
				for (j = 0; j < t -> w; j++) { /* scorrimento colonne */
					float x;
					x = get_normal_random(mean, std); /* calcolo gaussiana */
					set_val(t, i, j, c, x);
				}
			}
		}
		
		/* Aggiornamento delle statistiche */
		compute_stats(t);
		
	}
	else{
		printf("Errore ip_mat_init_random!!!\n");
		exit(1);
	}
}

/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat * ip_mat_copy(ip_mat * in) {
	
	if (in) {
		
		unsigned int i, j, c;
		
		/* Creazione della nuova ip_mat */
		ip_mat *out;
		out = ip_mat_create(in -> h, in -> w, in -> k, 0.0);
		
		for (i = 0; i < out -> h; i++) { /* scorrimento righe */
			for (j = 0; j < out -> w; j++) { /* scorrimento colonne */
				for (c = 0; c < out -> k; c++) { /* scorrimento canali */
					set_val(out, i, j, c, get_val(in, i, j, c));
				}
			}
		}
		
		/* Aggiornamento delle statistiche */
		compute_stats(out);
		
		return out;
	}
	else{
		printf("Errore ip_mat_copy!!!\n");
		exit(1);
	}
}

/* Restituisce una sotto-matrice */
ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end) {
	
	if (t && row_start < row_end && row_end <= t -> h && col_start < col_end && col_end <= t -> w) {
		
		unsigned int i, j, c;
		
		ip_mat *subset;
		
		/* Dimensioni della sotto-matrice*/
		unsigned int h, w;
		
		h = row_end-row_start;
		w = col_end-col_start;
		
		/* Creazione della sotto-matrice */
		subset = ip_mat_create(h, w, t -> k, 0.0);
		
		/* Loop */
		for (i = row_start; i < row_end; i++) { /* scorrimento righe */
			for (j = col_start; j < col_end; j++) { /* scorrimento colonne */
				for (c = 0; c < subset -> k; c++) { /* scorrimento canali */
					set_val(subset, i - row_start, j - col_start, c, get_val(t, i, j, c));
				}
			}
		}
		
		/* Aggiornamento delle statistiche */
		compute_stats(subset);
		
		return subset;		
	}
	else {
		printf("Errore ip_mat_subset!!!\n");
		exit(1);
	}
}

/* Concatena due ip_mat su una certa dimensione */
ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione){
	
	if (a && b) {
		
		unsigned int i, j, c;
		
		ip_mat* conc;
		
		/* valutazioni dimensione da concatenare */
		switch (dimensione) {
			case 0: { /* Caso h */
				if (a -> w == b -> w && a -> k == b -> k) {
					unsigned int h;
					h = a -> h + b -> h;
					conc = ip_mat_create(h, a -> w, a -> k, 0.0);
					for (i = 0; i < h; i++) {
						for (j = 0; j < a -> w; j++) {
							for (c = 0; c < a -> k; c++) {
								if (i < a -> h)
									set_val(conc, i, j, c, get_val(a, i, j, c));
								else
									set_val(conc, i, j, c, get_val(b, i - a -> h, j, c));
							}
						}
					}
				}
				else {
					printf("Errore ip_mat_concat!!!\n");
					exit(1);
				}
				break;
				
			}
			case 1: { /* Caso w */
				if (a -> h == b -> h && a -> k == b -> k) {
					unsigned int w;
					w = a -> w + b -> w;
					conc = ip_mat_create(a -> h, w, a -> k, 0.0);
					for (i = 0; i < a -> h; i++) {
						for (j = 0; j < w; j++) {
							for (c = 0; c < a -> k; c++){
								if (j < a -> w)
									set_val(conc, i, j, c, get_val(a, i, j, c));
								else
									set_val(conc, i, j ,c, get_val(b, i, j - a -> w, c));
							}
						}
					}
				}
				else {
					printf("Errore ip_mat_concat!!!\n");
					exit(1);
				}
				break;
			}
			case 2: { /* Caso k */
				if (a -> h == b -> h && a -> w == b -> w) {
					unsigned int k;
					k = a -> k + b -> k;
					conc = ip_mat_create(a -> h, a -> w, k, 0.0);
					for(i = 0; i < a -> h; i++) {
						for (j = 0; j < a -> w; j++) {
							for (c = 0; c < k; c++) {
								if (c < a -> k)
									set_val(conc, i, j, c, get_val(a, i, j, c));
								else
									set_val(conc, i, j, c, get_val(b, i, j, c - a -> k));
							}
						}
					}
				}
				else {
					printf("Errore ip_mat_concat!!!\n");
					exit(1);
				}
				break;
			}
			default: { /* Dimensione non valida */
				printf("Errore ip_mat_concat!!!\n");
				exit(1);
				break;
			}
		}
		
		/* Aggiornamento delle statistiche */
		compute_stats(conc);
		
		return conc;
	}
	else {
		printf("Errore ip_mat_concat!!!\n");
		exit(1);
	}
	
}

/**** PARTE 1: OPERAZIONI TRA IP_MAT ****/

/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output. */
ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b) {
	
	if (a && b && a -> h == b -> h && a -> w == b -> w && a -> k == b -> k) {
		
		unsigned int i, j, c;
		
		ip_mat *sum;
		sum = ip_mat_create(a -> h, a -> w, a -> k, 0.0);
		
		for (i = 0; i < sum -> h; i++) { /* scorrimento righe */
			for (j = 0; j < sum -> w; j++) { /* scorrimento colonne */
				for (c = 0; c < sum -> k; c++) { /* scorrimento canali */
					float temp;
					temp = get_val(a, i, j, c) + get_val(b, i, j, c);
					set_val(sum, i, j, c, temp);
				}
			}
		}
		
		/* Aggiornamento delle statische */
		compute_stats(sum);
		
		return sum;
		
	}
	else {
		printf("Errore ip_mat_sum!!!\n");
		exit(1);
	}
}

/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output */
ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b) {
	
	if	(a && b && a -> h == b -> h && a -> w == b -> w && a -> k == b -> k) {
		
		unsigned int i, j, c;
		ip_mat *sub;
		
		sub = ip_mat_create(a -> h, a -> w, a -> k, 0.0);
		
		for (i = 0; i < a -> h; i++){ /* scorrimento righe */
			for(j = 0; j < a -> w; j++) { /* scorrimento colonne */
				for (c = 0; c < a -> k; c++) { /* scorrimento canali */
					float temp;
					temp = get_val(a, i, j, c) - get_val(b, i, j, c);
					set_val(sub, i, j, c, temp);
				}
			}
		}
		
		/* Aggiornamento delle statistiche */
		compute_stats(sub);
		
		return sub;
	}
	
	else{
		printf("Errore ip_mat_sub!!!\n");
		exit(1);
	}
}

/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 * e si salva il risultato in un nuovo tensore in output. */
ip_mat * ip_mat_mul_scalar(ip_mat *a, float c) {
	
	if (a) {
		
		unsigned int i, j, l;
		
		ip_mat *mult;
		mult = ip_mat_create(a -> h, a -> w, a -> k, 0.0);
		
		for (i = 0; i < mult -> h; i++) { /* scorrimento righe */
			for (j = 0; j < mult -> w; j++) { /* scorrimento colonne */
				for (l = 0; l < mult -> k; l++) { /* scorrimento canali */
					float temp;
					temp = get_val(a, i, j, l) * c;
					set_val(mult, i, j, l, temp);
				}
			}
		}
		
		/* Aggiornamento delle statistiche */
		compute_stats(mult);
		
		return mult;
	}
	else{
		printf("Errore ip_mat_mul_scalar!!!\n");
		exit(1);
	}
}

/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output. */
ip_mat * ip_mat_add_scalar(ip_mat *a, float c) {
	
	if (a) {
		
		unsigned int i, j, l;
		
		ip_mat *add;
		add = ip_mat_create(a -> h, a -> w, a -> k, 0.0);
		
		for (i = 0; i < add -> h; i++) { /* scorrimento righe */
			for (j = 0; j < add -> w; j++) { /* scorrimento colonne */
				for (l = 0; l < add -> k; l++) { /* scorrimento canali */
					float temp;
					temp = get_val(a, i, j, l) + c;
					set_val(add, i, j, l, temp);
				}
			}
		}
		
		/* Aggiornamento delle statistiche */ 
		compute_stats(add);
		
		return add;
	}
	else{
		printf("Errore ip_mat_add_scalar!!!\n");
		exit(1);
	}
}


/* Restituisce la media di due matrici */
ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b) {
	
	if (a && b && a -> h == b -> h && a -> w == b -> w && a -> k == b -> k) {
		
		unsigned int i, j, c;
		
		ip_mat *mean;
		mean = ip_mat_create(a -> h, a -> w, a -> k, 0.0);
		
		for (i = 0; i < mean -> h; i++) { /* scorrimento righe */
			for (j = 0; j < mean -> w; j++) { /* scorrimento colonne */
				for (c = 0; c < mean -> k; c++) { /* scorrimento canali */
					float temp;
					temp = (get_val(a, i, j, c) + get_val(b, i, j, c)) / 2.0;
					set_val(mean, i, j, c, temp);
				}
			}
		}
		
		/* Aggiornamento delle statistiche */
		compute_stats(mean);
		
		return mean;
		
	}
	else {
		printf("Errore ip_mat_mean!!!\n");
		exit(1);
	}
}

/**** PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI ****/

/* Converte un'immagine RGB ad una immagine a scala di grigio */
ip_mat * ip_mat_to_gray_scale(ip_mat * in) {
	
	if (in) {
		
		unsigned int i, j, c;
		
		ip_mat * gray;
		gray = ip_mat_copy(in);
		
		for (i = 0; i < gray -> h; i++) { /* scorrimento righe */
			for (j = 0; j < gray -> w; j++) { /* scorrimento colonne */
				float acc = 0.;
				for (c = 0; c < gray -> k; c++) { /* scorrimento canali */
					acc += get_val(in, i, j, c);
				}
				acc /= gray -> k; /* media del canale */
				for (c = 0; c < gray -> k; c++) { /* scorrimento canali */
					set_val(gray, i, j, c, acc);
				}
			}
		}
		
		/* Aggiornamento delle statistiche */
		compute_stats(gray);
		
		return gray;
		
	}
	else{
		printf("Errore ip_mat_to_gray_scale!!!\n");
		exit(1);
	}
	
	
}

/* Effettua la fusione (combinazione convessa) di due immagini */
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha) {
	
	/* Controllo se le due ip_mat hanno la stessa dimensione */
	if (a && b && alpha >= 0 && alpha <= 1 && a -> h == b -> h && a -> w == b -> w && a -> k == b -> k) {
		
		/* ip_mat ausiarie */
		ip_mat *temp1, *temp2, *blend;
		
		temp1 = ip_mat_mul_scalar(a, alpha);
		temp2 = ip_mat_mul_scalar(b, (1 - alpha));
		blend = ip_mat_sum(temp1, temp2);
		
		/* Free delle ip_mat ausiliarie */
		ip_mat_free(temp1);
		ip_mat_free(temp2);
		
		/* Aggiornamento delle statistiche */
		compute_stats(blend);
		
		return blend;
		
	}
	else {
		printf("Errore ip_mat_blend!!!\n");
		exit(1);
	}
}

/* Operazione di brightening: aumenta la luminosità dell'immagine
 * aggiunge ad ogni pixel un certo valore */
ip_mat * ip_mat_brighten(ip_mat * a, float bright) {
	
	if (a) {
		ip_mat * br;
		
		br = ip_mat_add_scalar(a, bright);
		
		/* Aggiornamento delle statistiche */
		compute_stats(br);
		
		return br;
	}
	else{
		printf("Errore ip_mat_brighten!!!\n");
		exit(1);
	}
}

/* Operazione di corruzione con rumore gaussiano */
ip_mat * ip_mat_corrupt(ip_mat * a, float amount) {
	
	if(a){
		
		ip_mat *corrupt, *temp;
		float mean, std;
		
		mean = 0.;
		std  = amount / 2.;
		
		temp = ip_mat_create(a->h, a->w, a->k, 0.);
		
		ip_mat_init_random(temp, mean, std);
		corrupt = ip_mat_sum(a, temp);
		
		/* Deallocazione delle ip_mat ausiliarie */
		ip_mat_free(temp);
		
		/* Aggiornamento delle statistiche */
		compute_stats(corrupt);
		
		return corrupt;
	}
	else{
		printf("Errore ip_mat_corrupt!!!\n");
		exit(1);
	}
}

/**** PARTE 3: CONVLUZIONE E FILTRI ****/

/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a". 
 * Assumiamo che le dimensioni del filtro siano dispari */

ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f) { 
	
	if (a && f && a->k == f->k) {
		
		ip_mat * padded, * new;
		unsigned int pad_h, pad_w;
		unsigned int i, j, c;
		
		/* calcolo dei bordi */
		pad_h = (f->h - 1) / 2;
		pad_w = (f->w - 1) / 2;
		
		padded = ip_mat_padding(a, pad_h, pad_w);
		new = ip_mat_create(a->h, a->w, a->k, 0.);
		
		/* applicazione del filtro sulla ip_mat */
		for (c = 0; c < new->k; c++) { 
			for (i = 0; i < new->h; i++) {
				for (j = 0; j < new->w; j++) {
					unsigned int x, y;
					float acc;
					acc = 0;
					/* scorrimento del filtro */
					for (x = 0; x < f->h; x++) {
						for (y = 0; y < f->w; y++) {
							float tmp_pad, tmp_f;
							tmp_pad = get_val(padded, i+x, j+y, c);
							tmp_f = get_val (f, x, y, c);
							acc += tmp_f * tmp_pad;
						}
					}
					set_val(new, i, j, c, acc);
				}
			}
		}
		/* Deallocazione delle ip_mat ausiliarie */
		ip_mat_free(padded);
		
		/* Aggiornamento delle statistiche */
		compute_stats(new);
		
		return new;
	} 
	else {
		printf("Errore convolve!!!\n");
		exit(1);
	}   
}


/* Aggiunge un padding all'immagine */
ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w) {
	if (a) {
		
		unsigned int i, j, c;
		
		ip_mat *pad;
		
		pad = ip_mat_create(a -> h + 2 * pad_h, a -> w + 2 * pad_w, a -> k, 0.);
		
		/* copia dei valori di a */
		for (i = 0; i < a -> h; i++) {
			for (j = 0; j < a -> w; j++) {
				for (c = 0; c < a -> k; c++) {
					float temp;
					temp = get_val(a, i, j, c);
					set_val(pad, pad_h + i, pad_w + j, c, temp);
				}
			}
		}
		
		/* Aggiornamento delle statische */
		compute_stats(pad);
		
		return pad;
	}
	else {
		printf("Errore ip_mat_padding!!!\n");
		exit(1);
	}
}

/**** FILTRI ****/

/* Crea un filtro di sharpening */
ip_mat * create_sharpen_filter() {
	ip_mat * new;
	unsigned int dim, c;
	dim = 3;
	
	new = ip_mat_create(dim, dim, dim, -1.);
	
	for (c = 0; c < dim; c++) {
		set_val(new, 0, 0, c, 0.);
		set_val(new, 2, 0, c, 0.);
		set_val(new, 0, 2, c, 0.);
		set_val(new, 2, 2, c, 0.);
		set_val(new, 1, 1, c, 5.);
	}
	
	return new;
}

/* Crea un filtro per rilevare i bordi */
ip_mat * create_edge_filter() {
	ip_mat * new;
	unsigned int dim, c;
	dim = 3;
	
	new = ip_mat_create(dim, dim, dim, -1.);
	
	for (c = 0; c < dim; c++) {
		set_val(new, 1, 1, c, 8.);
	}
	
	return new;
}

/* Crea un filtro per aggiungere profondità*/
ip_mat * create_emboss_filter() {
	ip_mat * new;
	unsigned int dim, c;
	dim = 3;
	
	new = ip_mat_create(dim, dim, dim, 1.);
	
	for (c = 0; c < dim; c++) {
		set_val(new, 0, 0, c, -2.);
		set_val(new, 0, 1, c, -1.);
		set_val(new, 0, 2, c,  0.);
		set_val(new, 1, 0, c, -1.);
		set_val(new, 2, 0, c,  0.);
		set_val(new, 2, 2, c,  2.);
	}
	
	return new;
}

/* Crea un filtro medio per la rimozione del rumore */
ip_mat * create_average_filter(unsigned int w, unsigned int h, unsigned int k) {
	
	ip_mat * new;
	float c;
	c = 1. / (w * h);
	
	new = ip_mat_create(w, h, k, c);
	
	return new;
}

/* Crea un filtro gaussiano per la rimozione del rumore */

ip_mat * create_gaussian_filter(unsigned int h, unsigned int w, unsigned int k, float sigma) {
	unsigned int i, j, c;
	int cx, cy;
	float x, y, gauss = 0.;
	float sum = 0.;
	ip_mat *gaussian;
	
	gaussian = ip_mat_create(w, h, k, 0.0);
	
	/* centro kernel */
	cx = h / 2; 
	cy = w / 2;
	
	for(i = 0; i < gaussian -> h; i++){
		for(j = 0; j < gaussian -> w; j++) {
			x = (int)i - cx; 
			y = (int)j - cy;
			gauss = (1.0 /(2.0 * PI * (pow(sigma,2)))) * expf(-1. * ((pow(x,2) + pow(y,2)) / (2.0 * (pow(sigma,2))))); /* calcolo guassiana */
			for(c = 0; c < gaussian -> k; c++) {                
				set_val(gaussian,i,j,c,gauss);
			}
			sum += gauss;
		}
	}
	
	
	for(i = 0; i < gaussian -> h; i++) {
		for(j = 0; j < gaussian -> w; j++){
			for(c = 0; c < gaussian -> k; c++){
				float temp;
				temp = get_val(gaussian,i,j,c) / sum; /* normalizzazione dei valori */
				set_val(gaussian,i,j,c,temp);
			}
		}
	}
	
	return gaussian;
}

/**** PARTE 4: CLAMP E RESCALE ****/

/* Effettua una riscalatura dei dati tale che i valori siano in [0,new_max] */
void rescale(ip_mat * t, float new_max) {
	
	if (t && new_max >= 0) {
		
		unsigned int i, j, c;
		
		float min, max;
		
		compute_stats(t);
		
		for (c = 0; c < t -> k; c++) {
			min = t -> stat[c].min;
			max = t -> stat[c].max;
			for (i = 0; i < t -> h; i++) {
				for (j = 0; j < t -> w; j++) {
					float val;
					val = ((get_val(t, i, j, c) - min) / (max - min)) * new_max;
					set_val(t, i, j, c, val);
				}
			}
		}
		/* Aggiornamento delle statistiche */
		compute_stats(t);
	}
	else {
		printf("Errore rescale!!!\n");
		exit(1);
	}
}

/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.*/
void clamp(ip_mat * t, float low, float high){
	
	if (t && low <= high) {
		
		unsigned int i,j,c;
		
		for(i = 0; i < t -> h; i++){
			for(j = 0; j < t -> w; j++){
				for(c = 0; c < t -> k; c++){
					float sup;
					sup = get_val(t,i,j,c);
					
					if(sup < low)
						set_val(t,i,j,c,low);
					
					if(sup > high)
						set_val(t,i,j,c,high);
				}
			}
		}
		/* Aggiornamento delle statistiche */
		compute_stats(t);
		
	}
	else {
		printf("Errore clamp!!!\n");
		exit(1);
	}
}
