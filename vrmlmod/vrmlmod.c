#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "datastruct.h"
#include "util.h"

#define VTIME 10
#define ZTIME 15
#define RTIME 30
#define GTIME 5
#define GLOWTIME 2

/*Viewing preferences*/
#define ORB_RAD 0.1
#define LEVEL_HEIGHT 1.0
#define TREE_WIDTH 8.0
#define BRANCH_GROWTH 1.0
#define TILT_ANGLE 0.262

#define DNA 0
#define RNA 1
#define AMINO 2
#define SCALE 3

void world_setup(FILE* vrml_file, int flag_b) {
   fprintf(vrml_file, "#VRML V2.0 utf8\n");
   
   fprintf(vrml_file, "\n### WORLD SETUP ###\n");
   fprintf(vrml_file, "WorldInfo { title \"VRML model of Multiple Alignment\" }\n");
   fprintf(vrml_file, "NavigationInfo { avatarSize [0, 0, 0]}\n");
   if (!flag_b) {
      fprintf(vrml_file, "Background { skyColor [0 0 0] } #BLACK\n");
   } else {
      fprintf(vrml_file, "Background {\n");
      fprintf(vrml_file, "skyColor [0 0 0, 0 0 0, .734 .56 .56, .8 .36 .36, .98 .98 .66]\n");
	  fprintf(vrml_file, "skyAngle [0.7, 1.2, 1.3, 1.57]\n");
      fprintf(vrml_file, "groundColor [0 0 0, .54 .27 .07, .87 .71 .53, 1 1 1]\n");
      fprintf(vrml_file, "groundAngle [ 1.4, 1.53, 1.57]\n");
	  fprintf(vrml_file, "}\n");
   }
} 

void rot_setup(FILE* vrml_file, double flag_t) {
   fprintf(vrml_file, "\n### ROTATION ###\n");
   fprintf(vrml_file, "# Rotation Interpolator #\n"); 
   fprintf(vrml_file, "DEF ROTATOR OrientationInterpolator {\n");
   fprintf(vrml_file, "   key [0, .2, .3, .5, .7, .8, 1]\n");
   fprintf(vrml_file, "   keyValue [0 1 0 1.57, 0 1 0 1.97, 0 1 0 4.31, 0 1 0 4.71, 0 1 0 5.11, 0 1 0 7.45, 0 1 0 7.85]\n"); 
   fprintf(vrml_file, "}\n");
   fprintf(vrml_file, "# Rotation Timer #\n");
   fprintf(vrml_file, "DEF RTIMER TimeSensor { cycleInterval %f loop TRUE }\n", RTIME/flag_t);
}

void glow_setup(FILE* vrml_file) {
   fprintf(vrml_file, "\n### GLOWING ###\n");
   fprintf(vrml_file, "# Glow Interpolator #\n");   
   fprintf(vrml_file, "DEF GLOWER ColorInterpolator { key [0, .5, 1] keyValue [0 0 0, 0.5 0.5 0.5, 0 0 0]}\n"); 
   fprintf(vrml_file, "# Glow Timer #\n");
   fprintf(vrml_file, "DEF GLOWTIMER TimeSensor { cycleInterval %d loop TRUE }\n", GLOWTIME);
}

double max(double a, double b, double c) {
   if (a > b) {
      if (a > c) {
	     return a;
	  } else 
	     return c;
   } else if (b > c) {
      return b;
   } else return c;
}

void view_setup(nwk_node* root, FILE* vrml_file, double flag_z, int flag_v, int flag_x, double flag_t, int flag_s, int flag_f) {
   double height = flag_z*tree_height(root)*LEVEL_HEIGHT;
   double width = flag_z*TREE_WIDTH;
   double depth = (flag_f ? ORB_RAD*root->seq->size : ORB_RAD*flag_x);
   double dim1 = max(2.0*width, 2.0*depth, 0);
   double dim2 = max(2.0*width, 2.0*height, 2.0*depth);
   
   fprintf(vrml_file, "\n### VIEW POINTS ###\n");
   
   if (flag_s) {
      fprintf(vrml_file, "DEF VIEWSIDE Viewpoint {\n");
      fprintf(vrml_file, "   position 0 0 %f\n", dim1);
      fprintf(vrml_file, "   orientation 0 0 1 1.57\n");
      fprintf(vrml_file, "   fieldOfView 1.0\n");
      fprintf(vrml_file, "}\n");
   } 
   
   fprintf(vrml_file, "DEF VIEWFRONT Viewpoint {\n");
   fprintf(vrml_file, "   position 0 0 %f\n", TREE_WIDTH );
   fprintf(vrml_file, "}\n");

   if (!flag_s) {
      fprintf(vrml_file, "DEF VIEWSIDE Viewpoint {\n");
      fprintf(vrml_file, "   position 0 0 %f\n", dim1);
      fprintf(vrml_file, "   orientation 0 0 1 1.57\n");
      fprintf(vrml_file, "   fieldOfView 1.0\n");
      fprintf(vrml_file, "}\n");
   } 

   /* Zooming Setup */
   fprintf(vrml_file, "DEF ZTIMER TimeSensor { cycleInterval %d }\n", ZTIME);

   fprintf(vrml_file, "DEF ZOOM_FRONT PositionInterpolator {\n");
   fprintf(vrml_file, "   key [0, 1]\n");
   fprintf(vrml_file, "   keyValue [0 0 %f,   0 %f %f] }\n", TREE_WIDTH, height/3.0, dim2);  /* .125*tan(TILT_ANGLE)*dim2, dim2); */

   fprintf(vrml_file, "DEF ZOOM_SIDE PositionInterpolator {\n");
   fprintf(vrml_file, "   key [0, 1]\n");
   fprintf(vrml_file, "   keyValue [0 0 %f,   0 %f %f] }\n", TREE_WIDTH, height/3.0, dim1); 

   fprintf(vrml_file, "ROUTE ZTIMER.fraction_changed TO ZOOM_FRONT.set_fraction\n");
   fprintf(vrml_file, "ROUTE ZTIMER.fraction_changed TO ZOOM_SIDE.set_fraction\n");
   fprintf(vrml_file, "ROUTE ZOOM_FRONT.value_changed TO VIEWFRONT.position\n");
   fprintf(vrml_file, "ROUTE ZOOM_SIDE.value_changed TO VIEWSIDE.position\n");
    
   if (flag_v) { 
   
   fprintf(vrml_file, "# View Point Timer #\n");
   fprintf(vrml_file, "DEF VTIMER TimeSensor { cycleInterval %f loop TRUE }\n", VTIME / flag_t);
   fprintf(vrml_file, "# View Point Script #\n");
   fprintf(vrml_file, "DEF CYCLE Script {\n");
   fprintf(vrml_file, "   eventIn SFTime cycleTime\n");
   fprintf(vrml_file, "   eventOut SFInt32 output\n");
   fprintf(vrml_file, "   eventOut SFBool set_0\n");
   fprintf(vrml_file, "   eventOut SFBool set_1\n");
   fprintf(vrml_file, "   url \"javascript:\n");
   fprintf(vrml_file, "      function initialize() {\n");
   fprintf(vrml_file, "         output = 0;\n");
   fprintf(vrml_file, "      }\n");
   fprintf(vrml_file, "      function cycleTime(value, time) {\n");
   fprintf(vrml_file, "         if (output == 0) {\n");
   fprintf(vrml_file, "            output = 1;\n");
   fprintf(vrml_file, "            set_1 = TRUE;\n");
   fprintf(vrml_file, "         } else {\n");
   fprintf(vrml_file, "            output = 0;\n");
   fprintf(vrml_file, "            set_0 = TRUE;\n");
   fprintf(vrml_file, "         }\n");
   fprintf(vrml_file, "      }\n");
   fprintf(vrml_file, "   \"\n");
   fprintf(vrml_file, "}\n");
   } 
}

void stop_setup(FILE* vrml_file) {
   fprintf(vrml_file, "# Freeze Script #\n");
   fprintf(vrml_file, "DEF END Script {\n");
   fprintf(vrml_file, "   eventIn SFTime stop\n");
   fprintf(vrml_file, "   eventOut SFBool freeze\n");
   fprintf(vrml_file, "   eventOut SFInt32 textflag\n");   
   fprintf(vrml_file, "   url \"javascript:\n");
   fprintf(vrml_file, "      function initialize() {\n");
   fprintf(vrml_file, "         freeze = TRUE;\n");
   fprintf(vrml_file, "         textflag = -1;\n");
   fprintf(vrml_file, "      }\n");
   fprintf(vrml_file, "      function stop(value, time) {\n");
   fprintf(vrml_file, "         freeze = FALSE;\n");
   fprintf(vrml_file, "         textflag = 0;\n");
   fprintf(vrml_file, "      }\n");
   fprintf(vrml_file, "   \"\n");
   fprintf(vrml_file, "}\n");   
   
   fprintf(vrml_file, "\n### STOP ROTATION ###\n");
   
   fprintf(vrml_file, "# Stop Script #\n");
   fprintf(vrml_file, "DEF STOPPER Script {\n");
   fprintf(vrml_file, "   eventIn SFBool isActive\n");
   fprintf(vrml_file, "   eventOut SFBool stop\n");
   fprintf(vrml_file, "   url \"javascript:\n");
   fprintf(vrml_file, "      function initialize() {\n");
   fprintf(vrml_file, "         stop = TRUE;\n");
   fprintf(vrml_file, "      }\n");
   fprintf(vrml_file, "      function isActive(value, time) {\n");
   fprintf(vrml_file, "         if (value == TRUE) {\n");
   fprintf(vrml_file, "            stop = FALSE;\n");
   fprintf(vrml_file, "         }\n");
   fprintf(vrml_file, "      }\n");
   fprintf(vrml_file, "   \"\n");
   fprintf(vrml_file, "}\n");   
   
   fprintf(vrml_file, "# Touch Sphere #\n");
   fprintf(vrml_file, "Group {\n");
   fprintf(vrml_file, "   children [\n");
   fprintf(vrml_file, "      DEF STOP TouchSensor {}\n");
   fprintf(vrml_file, "      Shape {\n");
   fprintf(vrml_file, "         appearance Appearance {\n");
   fprintf(vrml_file, "            material Material {\n");
   fprintf(vrml_file, "               transparency 1.0\n");
   fprintf(vrml_file, "            }\n");
   fprintf(vrml_file, "         }\n");
   fprintf(vrml_file, "         geometry Sphere { radius .5 }\n");
   fprintf(vrml_file, "      }\n");
   fprintf(vrml_file, "   ]\n");
   fprintf(vrml_file, "}\n");
}

void drive_setup(nwk_node* n_node, FILE* vrml_file, double flag_t) {
   int i;
   int height = tree_height(n_node);
   fprintf(vrml_file, "\n### ANIMATION DRIVER ###\n");
   
   fprintf(vrml_file, "# Driver Script #\n");
   fprintf(vrml_file, "DEF DRIVE Script {\n");
   fprintf(vrml_file, "   eventIn SFTime cycleTime\n");
   for (i = 0; i < height+1; i++) {
      fprintf(vrml_file, "   eventOut SFTime go_%d\n", i);
   }
   fprintf(vrml_file, "   eventOut SFInt32 counter\n");
   fprintf(vrml_file, "   url \"javascript:\n");
   fprintf(vrml_file, "      function initialize() {\n");
   fprintf(vrml_file, "         counter = -1;\n");
   fprintf(vrml_file, "      }\n");
   fprintf(vrml_file, "      function cycleTime(value, time) {\n");
   for (i = 0; i < height+1; i++) {
      if (i == 0) { 
	     fprintf(vrml_file, "         if (counter == %d) {\n", i);
 	 } else {
	     fprintf(vrml_file, "         } else if (counter == %d) {\n", i);	 
	 }
	 fprintf(vrml_file, "            counter = %d;\n", i+1 );
	 fprintf(vrml_file, "            go_%d = value;\n", i); 
   }
   fprintf(vrml_file, "         } else {\n");
   fprintf(vrml_file, "            counter = counter + 1;\n");			
   fprintf(vrml_file, "         }\n");
   fprintf(vrml_file, "      }\n");
   fprintf(vrml_file, "   \"\n");
   fprintf(vrml_file, "}\n");   
   
   fprintf(vrml_file, "# Growth Timers #\n");
   fprintf(vrml_file, "DEF MASTERGROWER TimeSensor { cycleInterval %f loop TRUE }\n", GTIME/flag_t);
   for (i = 0; i < height; i++) {
      fprintf(vrml_file, "DEF GROWER%d TimeSensor { cycleInterval %f}\n", i, GTIME/flag_t);
   }
}

void dummy_branch(nwk_node* n_node, FILE* vrml_file, char* parent) {
   int i;
   char* new_parent = malloc((2+strlen(parent))*sizeof(char));
   /* to hold the string "X_\0" */

   for (i = 0; i < tg_size(n_node->children); i++) {
      fprintf(vrml_file, "   DEF BR_%s%d IndexedLineSet {\n", parent, i);
	  fprintf(vrml_file, "      color Color { color [.25 .25 .25, .25 .25 .25] }\n");
	  fprintf(vrml_file, "      coord DEF POINT_%s%d Coordinate { point [0 0 0, 0 0 0] }\n", parent, i);
	  fprintf(vrml_file, "      coordIndex [ 0 1 -1] }\n");
	  sprintf(new_parent, "%s%d_", parent, i);       
	  dummy_branch(n_node->children[i], vrml_file, new_parent);
   }
   free(new_parent);
}

char* seq_table(int type, int index) {
   if (type == DNA) {
      switch (index) {
	     case 0 : return "A";
		 case 1 : return "C";
		 case 2 : return "G";
		 case 3 : return "T";
		 default : 
		    printf("Invalid Character in Sequence\n");
	        return "";
	  }
   } else if (type == RNA) {
      switch (index) {
	     case 0 : return "A";
		 case 1 : return "C";
		 case 2 : return "G";
		 case 3 : return "U";
		 default : 
		    printf("Invalid Character in Sequence\n");
	        return "";
	  }
   } else if (type == AMINO) {
      switch (index) {
	     case 0  : return "A";
		 case 1  : return "C";
		 case 2  : return "D";
		 case 3  : return "E";
		 case 4  : return "F";
		 case 5  : return "G";
		 case 6  : return "H";
		 case 7  : return "I";
		 case 8  : return "K";
		 case 9  : return "L";
		 case 10 : return "M";
		 case 11 : return "N";
		 case 12 : return "P";
		 case 13 : return "Q";
		 case 14 : return "R";
		 case 15 : return "S";
 		 case 16 : return "T";
		 case 17 : return "V";
		 case 18 : return "W";
 		 case 19 : return "Y";
		 default : 
		    printf("Invalid Character in Sequence\n");
	        return "";
	  } 
   } else { /* type == scale */
      switch (index) {
	     case 0  : return "0";
		 case 1  : return "1";
		 case 2  : return "2";
		 case 3  : return "3";
		 case 4  : return "4";
		 case 5  : return "5";
		 case 6  : return "6";
		 case 7  : return "7";
		 case 8  : return "8";
		 case 9  : return "9";
	     case 10  : return "10";
		 case 11  : return "11";
		 case 12  : return "12";
		 case 13  : return "13";
		 case 14  : return "14";
		 case 15  : return "15";
		 case 16  : return "16";
		 case 17  : return "17";
		 case 18  : return "18";
		 case 19  : return "19";		 
		 default : 
		    printf("Invalid Character in Sequence\n");
	        return "";		 
	  }
   }
   return "";
}
		 
double* color_map(int type, int index) {
   double a;
   double b;
   double c; 
   double* ret = (double*) malloc(3*sizeof(double));
   if ((type == DNA) || (type == RNA)) {
      switch (index) {
	     case 0 : a = 0.5; b = 1.0; c = 0.2; break;
		 case 1 : a = 0.2; b = 0.5; c = 1.0; break;
		 case 2 : a = 1.0; b = 1.0; c = 0.2; break;
		 case 3 : 
	        if (type == DNA) {
			   a = 1.0; 
			   b = 0.5; 
			   c = 0.2;
			} else {
			   a = 1.0;
			   b = 0.2;
			   c = 0.2;
		   }
	  }
   } else if (type == AMINO) { /*mapping for amino acids*/
      switch (index % 3) {
	     case 0 : a = 0.25; break;
		 case 1 : a = 0.50; break;
		 case 2 : a = 0.75; break;
	  }
      switch (index % 4) {
	     case 0 : b = 0.2; break;
		 case 1 : b = 0.4; break;
		 case 2 : b = 0.6; break;
	     case 3 : b = 0.8; break;
	  }	  	  
      switch (index % 5) {
	     case 0 : c = 0.167; break;
		 case 1 : c = 0.333; break;
		 case 2 : c = 0.500; break;
	     case 3 : c = 0.666; break;
	     case 4 : c = 0.833; break;		 
	  }
   } else { /*mapping for scaled sequences*/
      switch (index) {
 	  	 case 0  : a = 0.10; b = 0.60; c = 0.1; break;
		 case 1  : a = 0.15; b = 0.60; c = 0.1; break;		 
		 case 2  : a = 0.20; b = 0.60; c = 0.1; break;
		 case 3  : a = 0.25; b = 0.60; c = 0.1; break;
		 case 4  : a = 0.30; b = 0.60; c = 0.1; break;
		 case 5  : a = 0.35; b = 0.60; c = 0.1; break;
 		 case 6  : a = 0.40; b = 0.60; c = 0.1; break;
		 case 7  : a = 0.45; b = 0.60; c = 0.1; break;
		 case 8  : a = 0.50; b = 0.60; c = 0.1; break;		 
		 case 9  : a = 0.55; b = 0.60; c = 0.1; break;
		 case 10 : a = 0.60; b = 0.60; c = 0.1; break;
		 case 11 : a = 0.60; b = 0.55; c = 0.1; break;
		 case 12 : a = 0.60; b = 0.50; c = 0.1; break;
         case 13 : a = 0.60; b = 0.45; c = 0.1; break;
		 case 14 : a = 0.60; b = 0.40; c = 0.1; break;
		 case 15 : a = 0.60; b = 0.35; c = 0.1; break;
         case 16 : a = 0.60; b = 0.30; c = 0.1; break;
		 case 17 : a = 0.60; b = 0.25; c = 0.1; break;
		 case 18 : a = 0.60; b = 0.20; c = 0.1; break;
		 case 19 : a = 0.60; b = 0.15; c = 0.1; break;
      }
   }
   ret[0] = a;
   ret[1] = b;
   ret[2] = c;
   return ret;
}

void dummy_orb(nwk_node* n_node, FILE* vrml_file, int flag_x, int flag_f) {
   int i, orbs = 0;
   int type = n_node->seq->type;
   double* color;
   if ((type == DNA) || (type == RNA)) {
      orbs = 4;
   } else { 
      orbs = 20;
   }
   if ((n_node->seq->size >  2 * flag_x) && !flag_f) {
      orbs = 20; /*twenty different shades for the amount of change */
	  type = SCALE;
   }
   fprintf(vrml_file, "# Dummy Orbs #\n");
   for (i = 0; i < orbs; i++) {
	  color = color_map(type,i);
	  fprintf(vrml_file, "   DEF ORB_%s Group {\n", seq_table(type, i));
	  fprintf(vrml_file, "      children [\n");
	  fprintf(vrml_file, "   Billboard {\n");
	  fprintf(vrml_file, "      children [\n");
	  fprintf(vrml_file, "         Transform {\n");
	  fprintf(vrml_file, "            translation 0 0 %f\n", ORB_RAD);
	  fprintf(vrml_file, "            children Text {\n");
	  fprintf(vrml_file, "               string \"%s\"\n", seq_table(type, i));
	  fprintf(vrml_file, "               fontStyle FontStyle { family \"TYPEWRITER\" size %f }}}]\n", ORB_RAD/2.0);
	  fprintf(vrml_file, "   }\n");
	  fprintf(vrml_file, "   Shape {\n");
	  fprintf(vrml_file, "      appearance Appearance {\n");
	  fprintf(vrml_file, "         material Material {\n");
  	  fprintf(vrml_file, "            diffuseColor %f %f %f } }\n", color[0], color[1], color[2]);
	  fprintf(vrml_file, "      geometry Sphere { radius %f } }\n", ORB_RAD);
	  fprintf(vrml_file, "   ] }\n");
   }
   fprintf(vrml_file, "# Dummy Glow Orbs #\n");
   for (i = 0; i < orbs; i++) {
      color = color_map(type,i);
      fprintf(vrml_file, "   DEF GORB_%s Group {\n", seq_table(type, i));
	  fprintf(vrml_file, "      children [\n");
	  fprintf(vrml_file, "   Billboard {\n");
	  fprintf(vrml_file, "      children [\n");
	  fprintf(vrml_file, "         Transform {\n");
	  fprintf(vrml_file, "            translation 0 0 %f\n", ORB_RAD);
	  fprintf(vrml_file, "            children [ Text {\n");
	  fprintf(vrml_file, "               string \"%s\"\n", seq_table(type, i));
	  fprintf(vrml_file, "               fontStyle FontStyle { family \"TYPEWRITER\" size %f }}]}]\n", ORB_RAD/2.0);
	  fprintf(vrml_file, "   }\n");
      fprintf(vrml_file, "   Shape {\n");
	  fprintf(vrml_file, "      appearance Appearance {\n");
	  fprintf(vrml_file, "         material DEF GLOW_%s Material {\n", seq_table(type, i));
  	  fprintf(vrml_file, "            diffuseColor %f %f %f } }\n", color[0], color[1], color[2]);
	  fprintf(vrml_file, "      geometry Sphere { radius %f } }\n", ORB_RAD);
	  fprintf(vrml_file, "   ] }\n");
   }
   fprintf(vrml_file, "# Dummy Transparent Orb #\n");
   fprintf(vrml_file, "   DEF TORB Shape {\n");
   fprintf(vrml_file, "      appearance Appearance {\n");
   fprintf(vrml_file, "         material DEF GLOW_TRANS Material {\n");
   fprintf(vrml_file, "            diffuseColor 0.0 0.0 0.0\n");
   fprintf(vrml_file, "            transparency 0.9 } }\n");
   fprintf(vrml_file, "      geometry Sphere { radius %f } }\n", ORB_RAD);
}


void proto_arc(nwk_node* n_node, FILE* vrml_file) {
   double i;
   
   fprintf(vrml_file, "PROTO ARC [\n");
   fprintf(vrml_file, "   field SFVec3f rad 1 1 1\n");
   fprintf(vrml_file, "]\n");
   fprintf(vrml_file, "{\n");
   fprintf(vrml_file, "   Transform {\n");
   fprintf(vrml_file, "      scale IS rad\n");
   fprintf(vrml_file, "      children [\n");
   fprintf(vrml_file, "      Shape {\n");
   fprintf(vrml_file, "      appearance Appearance {\n");
   fprintf(vrml_file, "         material Material {\n");
   fprintf(vrml_file, "            diffuseColor .25 .25 .25 \n");
   fprintf(vrml_file, "            emissiveColor .5 .5 .5 } }\n");
   fprintf(vrml_file, "      geometry IndexedLineSet {\n");
   fprintf(vrml_file, "         coord Coordinate {\n");
   fprintf(vrml_file, "         point [ ");
   for (i = 0; i <= 3.14; i += .098125) {
      if (i == 0) {
	     fprintf(vrml_file, "1 0 0");
	  } else {
	     fprintf(vrml_file, ", %f 0 %f", cos(i), sin(i));
      }
   }
   fprintf(vrml_file, "] }\n");   
   fprintf(vrml_file, "         coordIndex [ "); 
   for (i = 0; i < 32; i++) {
      fprintf(vrml_file, "%d %d -1 ", (int)i, (int)(i+1));
   }
   fprintf(vrml_file, "]\n"); 
   
   fprintf(vrml_file, "      }} \n");
   fprintf(vrml_file, "   ]}}\n");
}
void dummy_setup(nwk_node* n_node, FILE* vrml_file, int flag_x, int flag_f) {
   fprintf(vrml_file, "\n### DUMMY NODES ###\n");
   fprintf(vrml_file, "Switch {\n");
   fprintf(vrml_file, "   choice [\n");
   dummy_orb(n_node, vrml_file, flag_x, flag_f);
   fprintf(vrml_file, "# Dummy Branches #\n");
   dummy_branch(n_node, vrml_file, "0_");
   fprintf(vrml_file, "# Dummy SS_Connect #\n");
   fprintf(vrml_file, "   ]\n");
   fprintf(vrml_file, "}\n");   
   if (flag_f || (n_node->seq->size <=  2 * flag_x) ) {
      proto_arc(n_node, vrml_file);
   }
}

/* we know the width to the left and the right */
double coord_map(nwk_node* n_node, int index, int sibs, double width) {
   double sib_distance = width / sibs;
   double coord = -0.5*width + 0.5*sib_distance;
   int i;

   /*coord map not called from root*/
   for (i = 0; i <index; i++) {
      coord += sib_distance;
   }    
  /* printf("Coordinates: Index- %d Sibs- %d Width- %f Coord- %f \n", index, sibs, width, coord);*/
   return sibs!=1 ? coord : 0;
}

void interpol_seq(nwk_node* n_node, FILE* vrml_file, char* parent, double width, double height, double flag_z) {
   int i;
   char* new_parent = malloc((2+strlen(parent))*sizeof(char));
   /* to hold the string "X_\0" */

   if (n_node->parent) { /*if not root run for yourself*/
      fprintf(vrml_file, "DEF PI_%s%d PositionInterpolator {\n", parent, n_node->index);
      fprintf(vrml_file, "   key [0, 1]\n");
	  fprintf(vrml_file, "   keyValue [0 0 0,   0 %f %f] }\n", 
	     height, coord_map(n_node, n_node->index, tg_size(n_node->parent->children), width)); 
   } 
   for (i = 0; i < tg_size(n_node->children); i++) {
   	  sprintf(new_parent, "%s%d_", parent, n_node->index);      
      interpol_seq(n_node->children[i], vrml_file, new_parent,  n_node->parent ? width/tg_size(n_node->parent->children) : width, height*BRANCH_GROWTH, flag_z); /* mistake */
   }
   free(new_parent);
}

void interpol_br(nwk_node* n_node, FILE* vrml_file, char* parent, double width, double height, double flag_z) {
   int i;
   char* new_parent = malloc((2+strlen(parent))*sizeof(char));
   /* to hold the string "X_\0" */

   if (n_node->parent) { /*if not root run for yourself*/
      fprintf(vrml_file, "DEF CI_%s%d CoordinateInterpolator {\n", parent, n_node->index);
      fprintf(vrml_file, "   key [0, 1]\n");
	  fprintf(vrml_file, "   keyValue [0 0 0, 0 0 0, 0 0 0, 0 %f %f] }\n", 
	     height, coord_map(n_node, n_node->index, tg_size(n_node->parent->children), width)); 
   } 
   for (i = 0; i < tg_size(n_node->children); i++) {
   	  sprintf(new_parent, "%s%d_", parent, n_node->index);      
      interpol_br(n_node->children[i], vrml_file, new_parent, n_node->parent ? width/tg_size(n_node->parent->children) : width, height*BRANCH_GROWTH, flag_z); /* mistake 2? */
   }
   free(new_parent);
}
    
void interpol_setup(nwk_node* n_node, FILE* vrml_file, double flag_z) {
   fprintf(vrml_file, "\n### INTERPOLATORS ###\n");
   fprintf(vrml_file, "# Position Interpolators for Orbs #\n");
   interpol_seq(n_node, vrml_file, "\0", TREE_WIDTH*flag_z, LEVEL_HEIGHT*flag_z, flag_z);
   fprintf(vrml_file, "\n# Coordinate Interpolators for Branches #\n");   
   interpol_br(n_node, vrml_file, "\0", TREE_WIDTH*flag_z, LEVEL_HEIGHT*flag_z, flag_z);
}

/*Tree Structure Creation*/

void tree_node_add(nwk_node* n_node, FILE* vrml_file, char* parent, int flag_z, int flag_a) {
   int i, j, k;
   double index;
   char* new_parent = malloc((2+strlen(parent))*sizeof(char));
   nwk_node* root = n_node;
   tg_string_vec seqs = n_node->seq->seqs;
   tg_string_vec parent_seqs;
   
   fprintf(vrml_file, "DEF SEQ_%s%d Transform{\n", parent, n_node->index);
   fprintf(vrml_file, "   children [\n");
   for (i = 0; i < tg_size(n_node->children); i++) {
      sprintf(new_parent, "%s%d_", parent, n_node->index);   
      tree_node_add(n_node->children[i], vrml_file, new_parent, flag_z, flag_a);
   }
   
   if (n_node->parent) {
      if (flag_a) {
	     while (root->parent) { root = root->parent; }
         parent_seqs = root->seq->seqs;
      } else { parent_seqs = n_node->parent->seq->seqs; }
      for (i = 0; i < tg_size(seqs); i++) {
         for (j = 0; j < strlen(seqs[i]); j++) {
		    index = 2*ORB_RAD*(i*strlen(seqs[0]) + j); 
	        if ((seqs[i][j] == '-') || (seqs[i][j] == '.')) { /*deletion*/
			   fprintf(vrml_file, "Transform { translation %f 0 0 children USE TORB }\n", index);
			} else if (seqs[i][j] != parent_seqs[i][j]) {
			   fprintf(vrml_file, "Transform { translation %f 0 0 children USE GORB_%c }\n", index, seqs[i][j]);
			} else { 
			   fprintf(vrml_file, "Transform { translation %f 0 0 children USE ORB_%c }\n", index, seqs[i][j]);
			}
		 }
	  }
   } else { /* case for root*/
      for (i = 0; i < tg_size(seqs); i++) {
         for (j = 0; j < strlen(seqs[i]); j++) {
		    index = 2*ORB_RAD*(i*strlen(seqs[0]) + j); 
	        if ((seqs[i][j] == '-') || (seqs[i][j] == '.')) { /*deletion*/
			   fprintf(vrml_file, "Transform { translation %f 0 0 children USE TORB }\n", index);
            } else {
 		       fprintf(vrml_file, "Transform { translation %f 0 0 children USE ORB_%c }\n", index, seqs[i][j]);
		    }
		 }
	  }
   }
   for (k = 0; k < tg_size(n_node->children); k++) {  
      fprintf(vrml_file, "# Branch %d #\n", k);
      for (i = 0; i < tg_size(seqs); i++) {
         for (j = 0; j < strlen(seqs[i]); j++) {
		    index = 2*ORB_RAD*(i*strlen(seqs[0]) + j); 
 		    fprintf(vrml_file, "Transform { translation %f 0 0 children USE BR_%s%d_%d }\n",
		       index, parent, n_node->index, k);
		 }
	  }
   }
  
   fprintf(vrml_file, "DEF TEXT_SWITCH%s%d Switch{\n", parent, n_node->index);
   fprintf(vrml_file, "   choice [\n");
   fprintf(vrml_file, "Transform {\n");
   fprintf(vrml_file, "   rotation 0 0 1 1.57\n");
   if ((tg_size(n_node->children) % 2) == 0) {      
      fprintf(vrml_file, "   translation %f %f %f\n", ORB_RAD*n_node->seq->size, flag_z*LEVEL_HEIGHT/4.0, 0.0);
   } else { 
      fprintf(vrml_file, "   translation %f %f %f\n", ORB_RAD*n_node->seq->size, flag_z*LEVEL_HEIGHT/4.0, 2*ORB_RAD);
   }
   fprintf(vrml_file, "      children [\n");
   fprintf(vrml_file, "         Billboard {\n");
   fprintf(vrml_file, "            axisOfRotation 1 0 0\n");
   fprintf(vrml_file, "            children Shape {\n");
   fprintf(vrml_file, "               appearance Appearance { material Material { diffuseColor 1.0 1.0 1.0 emissiveColor 1.0 1.0 1.0} } \n");
   fprintf(vrml_file, "               geometry Text {\n");
   fprintf(vrml_file, "                  string \"%s\"\n", n_node->label);
   fprintf(vrml_file, "                  fontStyle FontStyle { family \"TYPEWRITER\" size %f style \"BOLD\" justify \"BEGIN\"}}}}]\n", 2.5 * ORB_RAD * flag_z);
   fprintf(vrml_file, "}]}\n");
   fprintf(vrml_file, "   ]\n");
   fprintf(vrml_file, "}\n");

}

int scale_color_map(int changes, int size_map, int seq_size) {
   int elmts_per_comp = (int) (seq_size/size_map);
   double ratio = (double) changes / (double) elmts_per_comp;
   if (ratio == 0) {
      return 0;
   } else if (ratio <= .01) {
      return 1;
   } else if (ratio <= .02) {
      return 2;
   } else if (ratio <= .04) {
      return 3;
   } else if (ratio <= .07) {
      return 4;
   } else if (ratio <= .10) {
      return 5;
   } else if (ratio <= .13) {
      return 6;
   } else if (ratio <= .16) {
      return 7;
   } else if (ratio <= .20) {
      return 8;
   } else if (ratio <= .24) {
      return 9;
   } else if (ratio <= .28) {
      return 10;
   } else if (ratio <= .32) {
      return 11;
   } else if (ratio <= .37) {
      return 12;
   } else if (ratio <= .43) {
      return 13;
   } else if (ratio <= .49) {
      return 14;
   } else if (ratio <= .55) {
      return 15;
   } else if (ratio <= .60) {
      return 16;
   } else if (ratio <= .70) {
      return 17;
   } else if (ratio <= .80) {
      return 18;
   } else {
      return 19;
   }
}
   
void scaled_tree_node_add(nwk_node* root, nwk_node* n_node, FILE* vrml_file, char* parent, int flag_x, int flag_z, int flag_a) {
   int i, j, k = 0, m = 0;
   int changes = 0;
   int changes_parent = 0;
   int own_color, parent_color;
   double index;
   char* new_parent = malloc((2+strlen(parent))*sizeof(char));
   tg_string_vec seqs = n_node->seq->seqs;
   tg_string_vec parent_seqs;
   
   fprintf(vrml_file, "DEF SEQ_%s%d Transform{\n", parent, n_node->index);
   fprintf(vrml_file, "   children [\n");
   for (i = 0; i < tg_size(n_node->children); i++) {
      sprintf(new_parent, "%s%d_", parent, n_node->index);   
      scaled_tree_node_add(root, n_node->children[i], vrml_file, new_parent, flag_x, flag_z, flag_a);
   }
   
   if (n_node->parent) {
      parent_seqs = (flag_a ? n_node->parent->seq->seqs : root->seq->seqs);
      for (i = 0; i < tg_size(seqs); i++) {
			   if (flag_x == m) {
			      break;
			   } 
	     for (j = 0; j < strlen(seqs[i]); j++) {
		    k++; /*keeps total objects we've examined so far*/
			if (seqs[i][j] != parent_seqs[i][j])  {
			   changes++;
		    }
			if (n_node->parent->seq->seqs[i][j] != parent_seqs[i][j])  {
			   changes_parent++;
		    }			
		    if ((k % ((int) root->seq->size/flag_x)) == 0) {
			   m++;
			   own_color = scale_color_map(changes, flag_x, root->seq->size);
			   parent_color = scale_color_map(changes_parent, flag_x, root->seq->size);
			   fprintf(vrml_file, "Transform { translation %f 0 0 children USE %s_%d }\n", 2*ORB_RAD*(m-1), 
			      ((own_color > (parent_color + 4)) ? "GORB" : "ORB"),
				   own_color);
			   if (flag_x == m) {
			      break;
			   } 
			   changes = 0;
			   changes_parent = 0;
			} 
		 }
	  }
   } else {  /* case for root*/
      for (i = 0; i < flag_x; i++) {
         index = 2.0*ORB_RAD*i; 
		 fprintf(vrml_file, "Transform { translation %f 0 0 children USE ORB_%d }\n", index, 0);
	  }
   }
   for (i = 0; i < tg_size(n_node->children); i++) {  
      fprintf(vrml_file, "# Branch %d #\n", i);
      for (j = 0; j < flag_x; j++) {
	     index = 2*ORB_RAD*j; 
		 fprintf(vrml_file, "Transform { translation %f 0 0 children USE BR_%s%d_%d }\n",
		       index, parent, n_node->index, i);
	  }
   }
   /*Name of Sequence*/
   fprintf(vrml_file, "DEF TEXT_SWITCH%s%d Switch{\n", parent, n_node->index);
   fprintf(vrml_file, "   choice [\n");
   fprintf(vrml_file, "Transform {\n");
   fprintf(vrml_file, "   rotation 0 0 1 1.57\n");
   fprintf(vrml_file, "   translation %f %f %f\n", ORB_RAD*flag_x, flag_z*LEVEL_HEIGHT/4.0, 2*ORB_RAD);
   fprintf(vrml_file, "      children [\n");
   fprintf(vrml_file, "         Billboard {\n");
   fprintf(vrml_file, "            axisOfRotation 1 0 0\n");
   fprintf(vrml_file, "            children Shape {\n");
   fprintf(vrml_file, "               appearance Appearance { material Material { diffuseColor 1.0 1.0 1.0 emissiveColor 1.0 1.0 1.0} } \n");
   fprintf(vrml_file, "               geometry Text {\n");
   fprintf(vrml_file, "                  string \"%s\"\n", n_node->label);
   fprintf(vrml_file, "                  fontStyle FontStyle { family \"TYPEWRITER\" size %f style \"BOLD\" justify \"BEGIN\"}}}}]\n", 2.0* ORB_RAD * flag_z);
   fprintf(vrml_file, "}]}\n");
   fprintf(vrml_file, "   ]\n");
   fprintf(vrml_file, "}\n");
}

void scaled_tree_setup(nwk_node* n_node, FILE* vrml_file, int flag_x, int flag_z, int flag_a) {
   fprintf(vrml_file, "\n### TREE SEQUENCE ###\n");
   fprintf(vrml_file, "Transform {\n"); /*new*/
   fprintf(vrml_file, "   rotation 1 0 0 %f\n", TILT_ANGLE); /*new*/
   fprintf(vrml_file, "   children\n"); /*new*/
   fprintf(vrml_file, "DEF TREE Transform {\n");
   fprintf(vrml_file, "   rotation 0 1 0 1.57\n");
   fprintf(vrml_file, "   children Transform {\n");
   fprintf(vrml_file, "      translation -%f 0 0\n", ORB_RAD*flag_x);
   fprintf(vrml_file, "      children [\n");
   scaled_tree_node_add(n_node, n_node, vrml_file, "\0", flag_x, flag_z, flag_a);
   fprintf(vrml_file, "      ]\n");
   fprintf(vrml_file, "   }\n");
   fprintf(vrml_file, "}\n");   
   fprintf(vrml_file, "}\n"); /*new*/
}

void rna_markup_setup(nwk_node* n_node, FILE* vrml_file) {
   tg_char_vec ss = n_node->seq->RNA_ss;
   char c;
   int i;
   double rad, index;
   ss_stack* s = new_ss_stack();
   ss_elmt* popped;

   if (n_node->seq->type == RNA) {
      for (i = 0; i < tg_size(ss)-1; i++) {
	     c = ss[i];
		 if (c == '<') {
		    push(s, new_ss_elmt(c, i));
		 } else if (c == '>') {
		    if (!isempty(s)) {
			   popped = pop(s);
			   rad = ORB_RAD*((double)i - (double)popped->index);
			   index = ORB_RAD*(2*popped->index) + rad;
			   fprintf(vrml_file, "Transform { translation %f 0 0 children ARC { rad %f %f %f} }\n", index, rad, rad, rad);
            } else {
			   printf("Mismatched SS_cons in Stockholm file at basepair: %d\n", i);
			   return;
			}
		 } else if (c != '.') {
		   printf("Illegal character '%c' in SS_cons of Stockholm file at basepair: %d\n", c, i);
		 }
      }
   }
}
   
   
void tree_setup(nwk_node* n_node, FILE* vrml_file, int flag_x, int flag_z, int flag_a, int flag_f) {
   if ((n_node->seq->size > 2*flag_x) && !flag_f) {
      scaled_tree_setup(n_node, vrml_file, flag_x, flag_z, flag_a);
	  return;
   }
   fprintf(vrml_file, "\n### TREE SEQUENCE ###\n");
   fprintf(vrml_file, "Transform {\n"); /*new*/
   fprintf(vrml_file, "   rotation 1 0 0 %f\n", TILT_ANGLE); /*new*/
   fprintf(vrml_file, "   children\n"); /*new*/
   fprintf(vrml_file, "DEF TREE Transform {\n");
   fprintf(vrml_file, "   rotation 0 1 0 1.57\n");
   fprintf(vrml_file, "   children Transform {\n");
   fprintf(vrml_file, "      translation -%f 0 0\n", ORB_RAD*(n_node->seq->size));
   fprintf(vrml_file, "      children [\n");
   tree_node_add(n_node, vrml_file, "\0", flag_z, flag_a);
   rna_markup_setup(n_node, vrml_file);
   fprintf(vrml_file, "      ]\n");
   fprintf(vrml_file, "   }\n");
   fprintf(vrml_file, "}\n");
   fprintf(vrml_file, "}\n"); /*new*/
}

void route_glow(nwk_node* n_node, FILE* vrml_file, int flag_x, int flag_f) {
   int i, orbs = 0;
   int type = n_node->seq->type;
   if ((type == DNA) || (type == RNA)) {
      orbs = 4;
   } else { 
      orbs = 20;
   }
   if ((n_node->seq->size > 2*flag_x) && !flag_f) {
      orbs = 20; /*ten different shades for the amount of change */
	  type = SCALE;
   }
   fprintf(vrml_file, "ROUTE GLOWTIMER.fraction_changed TO GLOWER.set_fraction\n");  
   for (i = 0; i < orbs; i++) {
      fprintf(vrml_file, "ROUTE GLOWER.value_changed TO GLOW_%s.emissiveColor\n", seq_table(type, i));
   }
   fprintf(vrml_file, "ROUTE GLOWER.value_changed TO GLOW_TRANS.emissiveColor\n");
}

void route_growers(nwk_node* n_node, FILE* vrml_file, char* parent) {
   int i;
   char* new_parent = malloc((2+strlen(parent))*sizeof(char));
      
   for (i = 0; i < tg_size(n_node->children); i++) {
      fprintf(vrml_file, "\n# Sequence Animation at Level %d #\n", n_node->height);
      fprintf(vrml_file, "ROUTE GROWER%d.fraction_changed TO PI_%s%d.set_fraction\n", n_node->height, parent, i);
      fprintf(vrml_file, "ROUTE GROWER%d.fraction_changed TO CI_%s%d.set_fraction\n", n_node->height, parent, i);
	  
	  fprintf(vrml_file, "ROUTE PI_%s%d.value_changed TO SEQ_%s%d.translation\n", parent, i, parent, i);
	  fprintf(vrml_file, "ROUTE CI_%s%d.value_changed TO POINT_%s%d.point\n", parent, i, parent, i);
	  
   }   
   for (i = 0; i < tg_size(n_node->children); i++) {
      sprintf(new_parent, "%s%d_", parent, i);
	  route_growers(n_node->children[i], vrml_file, new_parent);
   }
}

void route_animate(nwk_node* n_node, FILE* vrml_file, int flag_r, int flag_s, int flag_v) {
   int i;
   
   fprintf(vrml_file, "\n# Sequence Animation Driver #\n");
   for (i = 0; i < tree_height(n_node); i++) {
      fprintf(vrml_file, "ROUTE DRIVE.go_%d TO GROWER%d.startTime\n", i, i);   
   }
   fprintf(vrml_file, "ROUTE DRIVE.go_0 TO ZTIMER.startTime\n");
   route_growers(n_node, vrml_file, "0_");
   if (!flag_r) {
      fprintf(vrml_file, "\n# Freeze after Growth #\n");
      fprintf(vrml_file, "ROUTE DRIVE.go_%d TO END.stop\n", i); /* Freeze after Growth */
      fprintf(vrml_file, "ROUTE END.freeze TO RTIMER.loop\n");
	  if (flag_v) {
  	     fprintf(vrml_file, "ROUTE END.freeze TO VTIMER.loop\n");
	     if (flag_s) { 
	        fprintf(vrml_file, "ROUTE END.freeze TO VIEWFRONT.set_bind\n");
         } else {
	        fprintf(vrml_file, "ROUTE END.freeze TO VIEWSIDE.set_bind\n");
         }
	  }
   }
}

void route_text(nwk_node* n_node, FILE* vrml_file, char* parent) {
   int i;
   char* new_parent = malloc((2+strlen(parent))*sizeof(char));

   for (i = 0; i < tg_size(n_node->children); i++) {
   	  fprintf(vrml_file, "ROUTE END.textflag TO TEXT_SWITCH%s%d.whichChoice\n", parent, i);
   }   
   for (i = 0; i < tg_size(n_node->children); i++) {
      sprintf(new_parent, "%s%d_", parent, i);
	  route_text(n_node->children[i], vrml_file, new_parent);
   }
}

void route_setup(nwk_node* n_node, FILE* vrml_file, int flag_r, int flag_x, int flag_v, int flag_s, int flag_f) {

   fprintf(vrml_file, "\n### ROUTING ###\n");
   fprintf(vrml_file, "\n# Rotation Routing #\n");
   
   if (!flag_r) {
      fprintf(vrml_file, "ROUTE RTIMER.fraction_changed TO ROTATOR.set_fraction\n");
      fprintf(vrml_file, "ROUTE ROTATOR.value_changed TO TREE.set_rotation\n");
	  fprintf(vrml_file, "\n# Stop Routing #\n");
      fprintf(vrml_file, "ROUTE STOP.isActive TO STOPPER.isActive\n");
      fprintf(vrml_file, "ROUTE STOPPER.stop TO RTIMER.enabled\n");
   }
   
   
   if (flag_v) {
      fprintf(vrml_file, "\n# Viewpoint Routing #\n");
      fprintf(vrml_file, "ROUTE VTIMER.cycleTime TO CYCLE.cycleTime\n");
      fprintf(vrml_file, "ROUTE CYCLE.set_0 TO VIEWSIDE.set_bind\n");  
      fprintf(vrml_file, "ROUTE CYCLE.set_1 TO VIEWFRONT.set_bind\n");
   }
   
   fprintf(vrml_file, "\n# Glow Routing #\n");
   route_glow(n_node, vrml_file, flag_x, flag_f);
   
   fprintf(vrml_file, "\n# Sequence Animation #\n");
   fprintf(vrml_file, "ROUTE MASTERGROWER.cycleTime TO DRIVE.cycleTime\n");
   route_animate(n_node, vrml_file, flag_r, flag_s, flag_v);

   fprintf(vrml_file, "ROUTE END.textflag TO TEXT_SWITCH%d.whichChoice\n", 0);
   route_text(n_node, vrml_file, "0_");
}
   
int model(nwk_node* root, FILE* vrml_file, int flag_r, double flag_t, int flag_x, 
          double flag_z, int flag_a, int flag_v, int flag_b, int flag_s, int flag_f) {
   
   world_setup(vrml_file, flag_b);
   if (!flag_r) {
      rot_setup(vrml_file, flag_t);
   }
   if (flag_z == 1.0) {
      if (tree_height(root) > 4) {
	     flag_z *= tree_height(root)-3;
	  }
   }
   
   glow_setup(vrml_file);
   view_setup(root, vrml_file, flag_z, flag_v, flag_x, flag_t, flag_s, flag_f);
   stop_setup(vrml_file);
   drive_setup(root, vrml_file, flag_t);
   dummy_setup(root, vrml_file, flag_x, flag_f);
   interpol_setup(root, vrml_file, flag_z);
   tree_setup(root, vrml_file, flag_x, flag_z, flag_a, flag_f);
   route_setup(root, vrml_file, flag_r, flag_x, flag_v, flag_s, flag_f);
   return 0;
}

/* argv = vrmlmod [options] [filename] [output] */
/* Flags so far:
   -r            -> no rotation
   -a            -> in short alignment mode  - compare sequences to ancestor and not parent for glow effect
                 -> in scaled alignment mode - compare sequences to parent and not ancestor for color code and glow effect
   -b            -> gradient background enable
   -t [double x] -> x times the normal rate
   -v            -> viewpoint enable
   -x [int x]    -> x objects to represent scaled down seq. default - 25  
   -z [int x]    -> x times the space 
*/
int main(int argc, char** argv) {
   FILE* vrml_file;
   nwk_node* root;
   int flag_a = 0; /*ancestor comparison*/
   int flag_r = 0; /* rotation */
   int flag_v = 0; /* enable viewpoint switching */
   int flag_x = 50; /* downsize factor */
   int flag_b = 0; /* background */
   int flag_s = 0; /* side-view becomes default */ 
   int flag_f = 0; /* override scaling */
   double flag_z = 1.0; /* space factor */
   double flag_t = 1.0; /* time factor */
   
   int args = argc-3;
   int i = 1;

   while (args > 0) {
      if (!strncmp("-", argv[i], 1)) {
	     switch (argv[i][1]) { 
		    case 'r' :
			   flag_r = 1;
			   i++;
			   args--;
			   break;
			case 't' :
			   i++;
			   flag_t = atof(argv[i]);
			   if ((flag_t > 10) || (flag_t <= 0)) { 
			      printf("Time option argument must be a floating point between 0 and 10\n");
				  return 1;
			   }
			   i++;
			   args -= 2;
			   break;
			case 'x' :
   		       i++;
			   flag_x = atoi(argv[i]);
			   if ((flag_x > 500) || (flag_x < 10)) { 
			      printf("Scale option argument must be an integer between 10 and 100\n");
				  return 1;
			   }
               i++;
			   args -= 2;
			   break;
			case 'z' :
			   i++;
			   flag_z = atof(argv[i]);
			   if ((flag_z > 5) || (flag_z <= 0)) { 
			      printf("Space option argument must be a floating point between 0.0 and 5.0\n");
				  return 1;
			   }
			   i++;
			   args -= 2;
			   break;
		    case 'a' :
			   flag_a = 1;
			   i++;
			   args--;
			   break;  
		    case 'v' :
			   flag_v = 1;
			   i++;
			   args--;
			   break;  
		    case 'b' :
			   flag_b = 1;
			   i++;
			   args--;
			   break;
			case 's' :
			   flag_s = 1;
			   i++;
			   args--;
			   break;
		    case 'f' :
			   flag_f = 1;
			   i++;
			   args--;
			   break;			   
		    default :
			   printf("Unrecognized Option: %s\n", argv[i]);
			   return 1;
	     }
	  } else {
	     printf("Too many arguments to function\n");
		 return 1;
	  }
   }
   if (argc >= 3) {
      root = nwk_parse(argv[argc-2]);
	  vrml_file = fopen(argv[argc-1], "w");
      if (root) {
	     model(root, vrml_file, flag_r, flag_t, flag_x, flag_z, flag_a, flag_v, flag_b, flag_s, flag_f);   
         return 0;
	  }
	  return 1;
   } else {
      printf("Too few arguments to function\n");
	  return 1;
   }
}
