/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/File_process.h"


int File_exists(const char *FileName){

  FILE *file;
 
  if (file = fopen(FileName, "r")){
        fclose(file);
        return 1;
  }

  return 0;
}


double get_units(char *str){
  if(strcmp(str, "mm") == 0){
      return Umm;
  }
  else if(strcmp(str, "cm") == 0){
      return Ucm;
  }
  else if(strcmp(str, "m") == 0){
      return Um;
  }
  else{
     	printf("\nWarning: Unknown units %s\n", str);
	return 0.0;
  }
}

int isUnsignedInt(char *str){
  if(str == NULL) return 0;

  int i=0;
  int len = strlen(str);

  for (i=0; i<len; i++){
    if (!isdigit(str[i]) && str[i] != 'e' && str[i] != 'E') return 0;
  }

  return 1;
}

int isUnsignedFloat(char *str){
  if(str == NULL) return 0;

  int withDecimal=0, withExp=0, i=0;
  int len = strlen(str);

  for (i=0; i<len; i++){
    if (!isdigit(str[i])){
      if (str[i] == '.'){
	if(withDecimal) return 0;
	withDecimal =1;
      }
      else if (str[i] == 'e' || str[i] == 'E'){
	if(withExp) return 0;
	withExp =1;
      }
      else if (str[i] == '-'){
        if(i > 0){
	  if(str[i-1] != 'e' && str[i-1] != 'E') return 0;
	}
	else return 0;
      }
      else if (str[i] == '+'){
        if(i > 0){
	  if(str[i-1] != 'e' && str[i-1] != 'E') return 0;
	}
	else return 0;
      }
      else return 0;
    }
  }

  return 1;
}

int isFloat(char *str){
  if(str == NULL) return 0;

  int withDecimal=0, withExp=0, i=0;
  int len = strlen(str);

  for (i=0; i<len; i++){
    if (!isdigit(str[i])){
      if (str[i] == '.'){
	if(withDecimal) return 0;
	withDecimal =1;
      }
      else if (str[i] == '-'){
        if(i > 0){
	  if(str[i-1] != 'e' && str[i-1] != 'E') return 0;
	}
      } 
      else if (str[i] == 'e' || str[i] == 'E'){
	if(withExp) return 0;
	withExp =1;
      }
      else if (str[i] == '+'){
        if(i > 0){
	  if(str[i-1] != 'e' && str[i-1] != 'E') return 0;
	}
	else return 0;
      }
      else return 0;
    }
  }

  return 1;
}

int isBoolean(char *str){
  if(str == NULL) return 0;

  else if(strcmp(str, "1") == 0 || strcmp(str, "0") == 0) return 1;
  else if(strcmp(str, "true") == 0 || strcmp(str, "false") == 0) return 1;
  else if(strcmp(str, "True") == 0 || strcmp(str, "False") == 0) return 1;
  else if(strcmp(str, "TRUE") == 0 || strcmp(str, "FALSE") == 0) return 1;
  else if(strcmp(str, "yes") == 0 || strcmp(str, "no") == 0) return 1;
  else if(strcmp(str, "Yes") == 0 || strcmp(str, "No") == 0) return 1;
  else if(strcmp(str, "YES") == 0 || strcmp(str, "NO") == 0) return 1;

  else return 0;
}


unsigned int getBoolean(char *str){
  if(str == NULL) return 2;

  else if(strcmp(str, "1") == 0 || strcmp(str, "true") == 0 || strcmp(str, "True") == 0 || strcmp(str, "TRUE") == 0 || strcmp(str, "yes") == 0 || strcmp(str, "Yes") == 0 || strcmp(str, "YES") == 0) return 1;

  else if(strcmp(str, "0") == 0 || strcmp(str, "false") == 0 || strcmp(str, "False") == 0 || strcmp(str, "FALSE") == 0 || strcmp(str, "no") == 0 || strcmp(str, "No") == 0 || strcmp(str, "NO") == 0) return 0;

  else return 2;
}


void CreateDir(char *DirName){

  #if defined(_MSC_VER)
    _mkdir(DirName, 0755);
  #else
    mkdir(DirName, 0755);
  #endif

}


int CopyFile(char *to, char *from){

  FILE *fd_to, *fd_from;
  char buf[4096];
  ssize_t nread;

  fd_from = fopen(from, "rb");
  if(fd_from == NULL) return -1;

  fd_to = fopen(to, "wb");
  if(fd_to == NULL){
    fclose(fd_from);
    return -1;
  }

  while (nread = fread(buf, 1, 4096, fd_from), nread > 0){
    char *out_ptr = buf;
    ssize_t nwritten;

    do{
      nwritten = fwrite(out_ptr, 1, nread, fd_to);
      if (nwritten >= 0){
        nread -= nwritten;
        out_ptr += nwritten;
      }
      else{
	fclose(fd_to);
	fclose(fd_from);
	return -1;
      }
    } while(nread > 0);
  }

  fclose(fd_to);
  fclose(fd_from);

  return 0;
}



void str_replace(char *search , char *replace , char *subject){

    char  *p = NULL , *old = NULL , *new_subject = NULL ;
    int c = 0 , search_size;
     
    search_size = strlen(search);
     
    //Count how many occurences
    for(p = strstr(subject , search) ; p != NULL ; p = strstr(p + search_size , search))
    {
        c++;
    }
     
    //Final size
    c = ( strlen(replace) - search_size )*c + strlen(subject);
     
    //New subject with new size
    new_subject = malloc( c );
     
    //Set it to blank
    strcpy(new_subject , "");
     
    //The start position
    old = subject;
     
    for(p = strstr(subject , search) ; p != NULL ; p = strstr(p + search_size , search))
    {
        //move ahead and copy some text from original subject , from a certain position
        strncpy(new_subject + strlen(new_subject) , old , p - old);
         
        //move ahead and copy the replacement text
        strcpy(new_subject + strlen(new_subject) , replace);
         
        //The new start position after this search match
        old = p + search_size;
    }
     
    //Copy the part after the last search match
    strcpy(new_subject + strlen(new_subject) , old);

    //Copy results
    strcpy(subject, new_subject);
     
    //free(new_subject);
}
