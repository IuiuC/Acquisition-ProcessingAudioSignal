#include <cvirte.h>		
#include <userint.h>
#include <ansi_c.h>
#include <formatio.h>
#include<advanlys.h>
#include <utility.h>
#include"proiectAPD.h"
#define SAMPLE_RATE		0
#define NPOINTS			1

int waveInfo[2];
double sampleRate = 0.0;
int npoints = 0;
int pointspersecond = 0;

double *waveData = 0;
static int panelHandle;
int *HistogramArray=0;
double *axisArray=0;

int main (int argc, char *argv[])
{
	if (InitCVIRTE (0, argv, 0) == 0)
		return -1;	/* out of memory */
	if ((panelHandle = LoadPanel (0, "proiectAPD.uir", PANEL)) < 0)
		return -1;
	DisplayPanel (panelHandle);
	RunUserInterface ();
	DiscardPanel (panelHandle);
	return 0;
}

double ZeroPoints()
{
	int k=0;
	for(int i=1;i<npoints;++i)
		if(waveData[i-1]*waveData[i]<0)
			k++;
	return k;		
}
int CVICALLBACK load (int panel, int control, int event,
					  void *callbackData, int eventData1, int eventData2)

{
	double maxVal = 0.0;
	double minVal = 0.0;
	double medie=0.0;
	double mediana=0.0;
	double dispersia=0.0;
	int maxIndex = 0;
	int minIndex = 0;
	int interv=20;
	double NzeroPoints=0;

	switch (event)
	{
		case EVENT_COMMIT:
			FileToArray ("C:\\faculta\\APD\\waveInfo.txt", waveInfo, VAL_INTEGER, 2, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			sampleRate=waveInfo[SAMPLE_RATE];
			npoints=waveInfo[NPOINTS];
			pointspersecond=npoints/6;
			
			waveData = (double *) calloc(npoints, sizeof(double));
			axisArray = (double *) calloc(interv, sizeof(double));
			
			 FileToArray ("C:\\faculta\\APD\\waveData.txt", waveData, VAL_DOUBLE, npoints, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			
			PlotY(panel, PANEL_GRAPH, waveData, npoints, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
			HistogramArray=(int *) calloc(npoints, sizeof(int));
			MaxMin1D(waveData,npoints,&maxVal,&maxIndex,&minVal,&minIndex);
			Mean(waveData,npoints,&medie);
			Median(waveData,npoints,&mediana);
			StdDev(waveData,npoints,&medie,&dispersia);
			Histogram(waveData, npoints, minVal-1, maxVal+1, HistogramArray, axisArray, interv); 
			NzeroPoints=ZeroPoints();
			
			
			SetCtrlVal(panel,PANEL_MIN,minVal);
			SetCtrlVal(panel,PANEL_MAX,maxVal);
			SetCtrlVal(panel,PANEL_MEDIE,medie);
			SetCtrlVal(panel,PANEL_MEDIANA,mediana);
			SetCtrlVal(panel,PANEL_DISP,dispersia);
			SetCtrlVal(panel, PANEL_ZERO, NzeroPoints);
			
			
			PlotXY(panel, PANEL_GRAPH_3, axisArray, HistogramArray, interv, VAL_DOUBLE, VAL_INTEGER, VAL_VERTICAL_BAR, VAL_SOLID_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_GREEN);
			break;
	}
	return 0;
}


