/* DEtxt3.h  Header file for text handling routines    
   Denton S. Ebel  last updated: 17-Apr-96  */

	void breakline1 (char buffer[], char part1[], char part2[], char part3[], char part4[]);
	void remblank (char buffer[]);
	void trimblank (char buffer[]);
	int getsubstring (char buffer[], char substring[], int firstchar, int lastchar);
	int getToChar (char buffer[], char substring[], int matchchar, int removstr);
	void trimslashn (char buffer[]);

/* end of DEtxt3.h */
