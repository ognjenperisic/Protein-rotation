/* Rotation of the molecule (pdb file), so its longest dimension is aligned with one
of the axes (X, Y or Z) */

#include "stdafx.h"
#include "windows.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <direct.h>
#include <cmath>
#include <math.h>
#include <conio.h>

using namespace std;

class p_atom
{
   public:

     p_atom* Load(char* cFileName);
     p_atom* NextLink(){ return this->link;};
     p_atom* NextLink_forward(){ return this->link_forward;};

     int Save(char* cFileName, char* extension);
     int ChainLink(p_atom* connect){this->link=connect; return 0;};
     int ChainLink_forward(p_atom* connect){this->link_forward=connect; return 0;};
     int Rotation(char what);
     
     int MinMax();
     
     int Translation();
     float Distance_3D();
     int Distance();

     float x,y,z;
     string description1;
     string description2;
          
     static p_atom* last_atom;
     static p_atom* First;
     static p_atom* max_1;  // two points
     static p_atom* max_2;  // which distance is maximum
     static p_atom* max_3;  // point with the maximum distance frow line (max_1 - max_2)
     static float min_x;
     static float min_y;
     static float min_z;

  private:
    p_atom* link;
    p_atom* link_forward;
};

void transforms(string a, float &br, int x, int y)
{
    char numb[256]="123.56";    
    char i;

    for (i=x;i<x+y;i++)
      numb[i-x]=a[i];
    
    br=atof(numb);
}

void transforms(string a, int &br, int x, int y)
{
    char numb[256]="";    
    char i;

    for (i=x;i<x+y;i++)
      numb[i-x]=a[i];
    
    br=atoi(numb);
}

int FileName(char* file_name,char* extension,char* original_file,char CurrProt)
{    
    int i=0;  
    file_name[0]='\0';
    do {				    
        file_name[i]=original_file[i];
        i++;
    }
    while (original_file[i]!='.');                    
    file_name[i++]='_';
    file_name[i++]=CurrProt;
    file_name[i]='\0';
    strcat(file_name,extension);                    

    return 0;
}

float p_atom::min_x=0;
float p_atom::min_y=0;
float p_atom::min_z=0;
p_atom* p_atom::last_atom;
p_atom* p_atom::First=NULL;
p_atom* p_atom::max_1;
p_atom* p_atom::max_2;
p_atom* p_atom::max_3;



p_atom* p_atom::Load(char* cFileName)
{       
    p_atom* carbon;
    ifstream fin(cFileName);    
    char readln[256];
    string line;
    float xc,yc,zc;
    BOOL begg=true;

    while (fin.getline(readln,256)) {
        line=readln;

        if (line.substr(0,4)=="ATOM")
        {
            if (line.substr(26,1)==" ") {                
                if (line.substr(13,2)!="HH") {
                    transforms(line,xc,27,12);
                    transforms(line,yc,39,8);
                    transforms(line,zc,47,8);
                    carbon=new p_atom;
                    carbon->ChainLink(First);
                    carbon->x=xc;
                    carbon->y=yc;
                    carbon->z=zc;
                    carbon->description1=line.substr(0,30);
                    carbon->description2=line.substr(55,25);
                    First=carbon;
                    
                    if (!begg)
                        carbon->NextLink()->ChainLink_forward(carbon);
                    if(begg) {
                        p_atom::last_atom=carbon;
                        begg=false;
                    }                    
                    First=carbon;               
                }
            } 
        }
        
        if (line.substr(0,6)=="HETATM")
        {

            if (line.substr(26,1)==" ") {                    
                
                if (line.substr(13,2)!="HH"){

            if (line.substr(26,1)==" ")
            {
                if (line.substr(13,2)!="HH")
                {

                    transforms(line,xc,27,12);
                    transforms(line,yc,39,8);
                    transforms(line,zc,47,8);
                    carbon=new p_atom;
                    carbon->ChainLink(First);
                    carbon->x=xc;
                    carbon->y=yc;
                    carbon->z=zc;
                    carbon->description1=line.substr(0,30);
                    carbon->description2=line.substr(55,25);
                    First=carbon;
                    
                    if (!begg)
                        carbon->NextLink()->ChainLink_forward(carbon);
                    if(begg) {
                        p_atom::last_atom=carbon;
                        begg=false;
                    }                    
                    First=carbon;               
                }
            } 
        }
    }   
    fin.close();
    carbon->ChainLink_forward(NULL);
    return carbon;
}

int p_atom::Save(char *cFileName, char* extension)
{
    p_atom* temp_atom;
    char file_name_temp[256];
    char file_name[256];

    file_name_temp[0]='\0';
    //strcpy(file_name_temp,dir1);
    //strcat(file_name_temp,"\\");
    strcat(file_name_temp,cFileName);
    
    FileName(file_name,extension,file_name_temp,'\0');
    ofstream fout(file_name);
    
    temp_atom=p_atom::last_atom;                   
    while(temp_atom) {
        fout<<temp_atom->description1;
        fout.width(8);
        fout.setf(ios::fixed);
        fout.precision(3);
        fout.setf(ios::right,ios::adjustfield);
        fout<<temp_atom->x;
        fout.width(8);
        fout.precision(3);
        fout.setf(ios::right,ios::adjustfield);
        fout<<temp_atom->y;
        fout.width(8);
        fout.precision(3);
        fout.setf(ios::right,ios::adjustfield);
        
        fout<<temp_atom->z;
        fout<<" "<<temp_atom->description2<<endl;
        
        temp_atom=temp_atom->NextLink_forward();
    }
    fout.close();
    return 0;
}

float p_atom::Distance_3D()
{
    p_atom* p0=this;

    float d_up;
    float d_down;
    float br_1_up,br_2_up,br_3_up;
    float br_1_down,br_2_down,br_3_down;

    /*             (y2       -        y1)        *       (z1         -    z0)           */
    br_1_up=((p_atom::max_2->y  -  p_atom::max_1->y)*(p_atom::max_1->z  -  p0->z)-
    /*             (z2       -        z1)        *       (y1         -    y0)           */
          (p_atom::max_2->z  -  p_atom::max_1->z)*(p_atom::max_1->y  -  p0->y));

    /*             (x2       -        x1)        *       (z1         -    z0)           */
    br_2_up=((p_atom::max_2->x  -  p_atom::max_1->x)*(p_atom::max_1->z  -  p0->z)-
    /*             (z2       -        z1)        *       (x1         -    x0)           */
          (p_atom::max_2->z  -  p_atom::max_1->z)*(p_atom::max_1->x  -  p0->x));

    /*             (x2       -        x1)        *       (y1         -     y0)           */
    br_3_up=((p_atom::max_2->x  -  p_atom::max_1->x)*(p_atom::max_1->y  -  p0->y)-
    /*             (y2       -        y1)        *       (x1         -    x0)           */
          (p_atom::max_2->y  -  p_atom::max_1->y)*(p_atom::max_1->x  -  p0->x));

    d_up=sqrt(br_1_up*br_1_up + br_2_up*br_2_up + br_3_up*br_3_up);

    br_1_down=((p_atom::max_2->x  -  p_atom::max_1->x));
    br_2_down=((p_atom::max_2->y  -  p_atom::max_1->y));
    br_3_down=((p_atom::max_2->z  -  p_atom::max_1->z));

    d_down=sqrt(br_1_down*br_1_down + br_2_down*br_2_down + br_3_down*br_3_down);

    if (d_down!=0)
        return d_up/d_down;
    else
        return 0;    
}

int p_atom::Distance()
{
    p_atom* temp_atom;
    float max_distance=0;
    float temp_distance;


    temp_atom=p_atom::last_atom;                   

    while(temp_atom) {
        temp_distance=temp_atom->Distance_3D();
        if (temp_distance>=max_distance) {
            p_atom::max_3=temp_atom;
            max_distance=temp_distance;
        }
        temp_atom=temp_atom->NextLink_forward();
    }
    return 0;
}

int p_atom::MinMax()
{
   p_atom* first_loop;
   p_atom* second_loop;
   float dist=0,max_dist=0;
   float dist1=0,dist2=0;
<<<<<<< HEAD
   

=======
>>>>>>> c7e77fa7c16a8e851af244176edeb731c5cbb669

   first_loop=p_atom::last_atom;
   
   while(first_loop) {                         
       
       second_loop=first_loop->NextLink_forward();

       while(second_loop) {
           dist=sqrt(pow(first_loop->x-second_loop->x,2)+pow(first_loop->y-second_loop->y,2)+pow(first_loop->z-second_loop->z,2));
           if (dist>max_dist) {
               dist1=sqrt(first_loop->x*first_loop->x+first_loop->y*first_loop->y+first_loop->z*first_loop->z);
               dist2=sqrt(second_loop->x*second_loop->x+second_loop->y*second_loop->y+second_loop->z*second_loop->z);
               max_dist=dist;

               if (dist2<=dist1){
                   p_atom::max_1=first_loop;
                   p_atom::max_2=second_loop;
               }
               else {
                   p_atom::max_2=first_loop;
                   p_atom::max_1=second_loop;
               }

           }
           
           second_loop=second_loop->NextLink_forward();
       }
       first_loop=first_loop->NextLink_forward();
   }

   
   return 0;
}


/*
int p_atom::MinMax_s()
{    
   p_atom* first_loop;
   p_atom* second_loop;
 
   float dist=0,max_dist=0;
   float dist1=0,dist2=0;
   first_loop=p_atom::last_atom;
   long int br=0;
   
   while(first_loop)
   {                         
       second_loop=p_atom::last_atom;
       //second_loop=first_loop->NextLink_forward();
       while(second_loop)
       {
           dist=sqrt(pow(first_loop->x-second_loop->x,2)+pow(first_loop->y-second_loop->y,2)+pow(first_loop->z-second_loop->z,2));
           if (dist>max_dist)
           {
               dist1=sqrt(first_loop->x*first_loop->x+first_loop->y*first_loop->y+first_loop->z*first_loop->z);
               dist2=sqrt(second_loop->x*second_loop->x+second_loop->y*second_loop->y+second_loop->z*second_loop->z);
               max_dist=dist;
               if (dist2<=dist1)
               {
                   p_atom::max_1=first_loop;
                   p_atom::max_2=second_loop;
               }
               else
               {
                   p_atom::max_2=first_loop;
                   p_atom::max_1=second_loop;
               }
           }
<<<<<<< HEAD
           br++;

=======
>>>>>>> c7e77fa7c16a8e851af244176edeb731c5cbb669
           second_loop=second_loop->NextLink_forward();
       }
       first_loop=first_loop->NextLink_forward();
   }

   cout<<"broj iteracija :"<<br<<endl;
   return 0;
}
*/


int p_atom::Translation()
{

    // translation
    // what == 0  XY plane
    // what == 1  XZ plane
    // what == 2  YZ plane

    p_atom* first_loop;    
    
    float min_x=p_atom::max_2->x;    
    float min_y=p_atom::max_2->y;    
    float min_z=p_atom::max_2->z;

   


    first_loop=p_atom::last_atom;
    while(first_loop) {
        first_loop->x=first_loop->x-min_x;
        first_loop->y=first_loop->y-min_y;
        first_loop->z=first_loop->z-min_z;
        
        first_loop=first_loop->NextLink_forward();
    }

    return 0;
}



int p_atom::Rotation(char what)
{

    // rotation
    // what == 0  XY plane
    // what == 1  XZ plane
    // what ==2   YZ plane
    p_atom* loop;
    loop=p_atom::last_atom;

    float max_x=p_atom::max_1->x,  min_x=p_atom::max_2->x;
    float max_y=p_atom::max_1->y,  min_y=p_atom::max_2->y;
    float max_z=p_atom::max_1->z,  min_z=p_atom::max_2->z;
    float x_temp,y_temp,z_temp;
    double sin_a,cos_a;
    double sin_b,cos_b;
    double R;
    float mp;
    
    if (!what) {
        mp=(-1)*((max_x)/fabs(max_x))*((max_y)/fabs(max_y));

        R=sqrt(pow(max_x-min_x,2)+pow(max_y-min_y,2)                     );
        sin_a=mp*abs(max_y-min_y)/R;  cos_a=(max_x-min_x)/R;  //rotation in XY plane
    }

    if (what==1) {
        mp=-1*((max_x)/fabs(max_x))*((max_z)/fabs(max_z));
        R=sqrt(pow(max_x-min_x,2)+                     pow(max_z-min_z,2));
        sin_b=mp*abs(max_z-min_z)/R;  cos_b=(max_x-min_x)/R;  //rotation in XZ plane
    }

    if (what==2){
        max_x=p_atom::max_3->x,  min_x=p_atom::max_3->x;
        max_y=p_atom::max_3->y,  min_y=0;
        max_z=p_atom::max_3->z,  min_z=0;

        mp=1*((max_y)/fabs(max_y))*((max_z)/fabs(max_z));
        R=                          sqrt(pow(max_y-min_y,2)+ pow(max_z-min_z,2));
        sin_b=mp*abs(max_y-min_y)/R;  cos_b=(max_z-min_z)/R;  //rotation in XZ plane
    }

    while(loop)
    {     
        
        
        x_temp=loop->x;
        y_temp=loop->y;
        z_temp=loop->z;
        
        if(!what) {
            loop->x=x_temp*cos_a-y_temp*sin_a;
            loop->y=x_temp*sin_a+y_temp*cos_a;
        }

        if(what==1) {
            loop->x=x_temp*cos_b-z_temp*sin_b;
            loop->z=x_temp*sin_b+z_temp*cos_b;
        }

        if(what==2) {
            loop->y=y_temp*cos_b-z_temp*sin_b;
            loop->z=y_temp*sin_b+z_temp*cos_b;
        }


        loop=loop->NextLink_forward();
    }
    return 0;
}


int main(int argc, char* argv[])
{
   
    
   WIN32_FIND_DATA FileData; 
   HANDLE hSearch;    
   BOOL fFinished = FALSE; 

   p_atom* carbon;
           
   if (argc < 2){
       printf("PDB file required.");
       exit(1);
   }
  
   hSearch = FindFirstFile(argv[1], &FileData); 
   if (hSearch == INVALID_HANDLE_VALUE){ 
       printf("No files found.");       
   } 
   else {	
       while (!fFinished) {  
           if (!fFinished) {
               printf("%s\n",FileData.cFileName);
               
               // START
               carbon=new p_atom;
               carbon=carbon->Load(FileData.cFileName);                    
               carbon->MinMax();
               carbon->Translation();

               carbon->Save(FileData.cFileName, "transl.pdb");

               carbon->Rotation(0);                    
               carbon->Rotation(1);

                    
               carbon->Save(FileData.cFileName, "test_pdb.trn");	

               carbon->Distance();
               carbon->Rotation(2);
                                                            
               
               carbon->Save(FileData.cFileName, "rot.pdb");				
               
               p_atom* pDel = carbon;
               while (pDel != NULL) {
                   
                   carbon = carbon->NextLink();
                   delete pDel;
                   pDel = NULL;
                   
                   pDel = carbon;
               }
               
               // FINISH
					
               fFinished=FALSE;
					 
               if (!FindNextFile(hSearch, &FileData)) {
                   if (GetLastError() == ERROR_NO_MORE_FILES) {                             
                       fFinished = TRUE; 
                       
                   } 
                   else { 
                       printf("Couldn't find next file."); 
                       return 0;
                   } 
               }
           }
       }
   }
	

   return 0;
}
