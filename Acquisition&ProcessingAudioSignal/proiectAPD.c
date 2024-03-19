#include <cvirte.h>		
#include <userint.h>
#include <ansi_c.h>
#include <formatio.h>
#include<advanlys.h>
#include <utility.h>
#include"proiectAPD.h"
#include <analysis.h>
#include "toolbox.h"
#define SAMPLE_RATE		0
#define NPOINTS			1

int secunda=1;
int N;
int waveInfo[2];
double sampleRate = 0.0;
int npoints = 0;
int start = 0;
int fs=44100;
double *waveData = 0;
double *waveDataD = 0;

double *waveCopy = 0;
double* waveCop=0;
static int panelHandle;
static int panelHandleFreq;

double* filteredData=0;
double* FilteredDataa=0;
int *HistogramArray=0;
double *axisArray=0;
double*deriv;
int startc=0, stopc=0;


WindowConst winConst;
WindowConst winConst1;
WindowConst winConstF;

	double *convertedSpectrum;
	double *convertedSpectrumF;
	
	double *autoSpectrum;
	double *autoSpectrumF;
	
	double df=0.0; 
	double freqPeak=0.0; 
	double powerPeak=0.0;
	
	  //voltage semnal

int main (int argc, char *argv[])
{
	
	int error = 0;
	
	/* initialize and load resources */
	nullChk (InitCVIRTE (0, argv, 0));
	errChk (panelHandle = LoadPanel (0, "proiectAPD.uir", PANEL));
	errChk (panelHandleFreq = LoadPanel (0, "proiectAPD.uir", PANEL_FREQ));
	 
	/* display the panel and run the user interface */
	errChk (DisplayPanel (panelHandle));
	 errChk (RunUserInterface ());

	Error:
	/* clean up */
	if (panelHandle > 0)
		DiscardPanel (panelHandle);
	return 0;
}


int CVICALLBACK OnSwitchPanelCB (int panel, int control, int event,
								 void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			if(panel == panelHandle)
			{
				SetCtrlVal(panelHandleFreq, PANEL_FREQ_BINARYSWITCH, 1);
				DisplayPanel(panelHandleFreq);
				HidePanel(panel);
			}
			else
			{
				SetCtrlVal(panelHandle,PANEL_BINARYSWITCH, 0);
				DisplayPanel(panelHandle);
				HidePanel(panel);
			}
			break;
	}
	return 0;
}
int CVICALLBACK On_Panel (int panel, int event, void *callbackData,
							   int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_GOT_FOCUS:

			break;
		case EVENT_LOST_FOCUS:

			break;
		case EVENT_CLOSE:
			QuitUserInterface (0);
			break;
	}
	return 0;
}


int CVICALLBACK On_Panel_Freq (int panel, int event, void *callbackData,
							   int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_GOT_FOCUS:
			
			break;
		case EVENT_LOST_FOCUS:

			break;
		case EVENT_CLOSE:
			QuitUserInterface (0); 
			break;
	}
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
	double asimetria;
	double aplatizare;



	
	switch (event)
	{
		case EVENT_COMMIT:
			FileToArray ("C:\\faculta\\APD_P\\waveInfo.txt", waveInfo, VAL_INTEGER, 2, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			sampleRate=waveInfo[SAMPLE_RATE];
			npoints=waveInfo[NPOINTS];
			//pointspersecond=npoints/6;
			
			waveData = (double *) calloc(npoints, sizeof(double));
			
			
			axisArray = (double *) calloc(interv, sizeof(double));
			
			 FileToArray ("C:\\faculta\\APD_P\\waveData.txt", waveData, VAL_DOUBLE, npoints, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			
			PlotY(panel, PANEL_GRAPH, waveData, npoints, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
			HistogramArray=(int *) calloc(npoints, sizeof(int));
			MaxMin1D(waveData,npoints,&maxVal,&maxIndex,&minVal,&minIndex);
			Mean(waveData,npoints,&medie);
			Median(waveData,npoints,&mediana);
			StdDev(waveData,npoints,&medie,&dispersia);
			Histogram(waveData, npoints, minVal-1, maxVal+1, HistogramArray, axisArray, interv); 
			NzeroPoints=ZeroPoints();
			//asimetria=Skewness();
			Moment(waveData, 256, 3, &asimetria);
			//aplatizare=Kurtosis();
			Moment(waveData, 256, 4, &aplatizare);
		
			
				
			SetCtrlVal(panel,PANEL_MIN,minVal);
			SetCtrlVal(panel,PANEL_MAX,maxVal);
			SetCtrlVal(panel,PANEL_MININDEX,minIndex);
			SetCtrlVal(panel,PANEL_MAXINDEX,maxIndex);
			SetCtrlVal(panel,PANEL_MEDIE,medie); 
			SetCtrlVal(panel,PANEL_MEDIANA,mediana);
			SetCtrlVal(panel,PANEL_DISP,dispersia);
			SetCtrlVal(panel, PANEL_ZERO, NzeroPoints);
			SetCtrlVal(panel,PANEL_SK, asimetria);
			SetCtrlVal(panel, PANEL_KTS, aplatizare);
			PlotXY(panel, PANEL_GRAPH_3, axisArray, HistogramArray, interv, VAL_DOUBLE, VAL_INTEGER, VAL_VERTICAL_BAR, VAL_SOLID_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);

			
			break;
	}
	return 0;
}

int CVICALLBACK On_Next (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	int *axaX=(int*)calloc(fs, sizeof(int));
	
	switch (event)
	{
		case EVENT_COMMIT:
			
		{
		if(start>=0 && start<=fs*5-1)
		{
			
			start+=fs;
			for(int i=0;i<fs;++i)
				{
					axaX[i]=start+i;
				}
			DeleteGraphPlot(panel, PANEL_GRAPH, -1, VAL_IMMEDIATE_DRAW);
			PlotXY(panel, PANEL_GRAPH, axaX, waveData+start, fs, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
			SetCtrlVal(panel,PANEL_START, start/fs);
			SetCtrlVal(panel,PANEL_STOP1, (start+fs)/fs);
			
			GetCtrlVal(panel,PANEL_STOP1, &stopc);
			GetCtrlVal(panel,PANEL_START, &startc);
			
		}
			
			break;
		}
	}
	return 0;
}

 

int CVICALLBACK On_Prev (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	int *axaX=(int*)calloc(fs, sizeof(int));
	
	switch (event)
	{
		case EVENT_COMMIT:
		{	
		
		//if(start>=fs && start<=fs*6-1)
		if(start>=fs)
		{
			start-=fs;
			for(int i=0;i<fs;++i)
			{
				axaX[i]=start+i;
			}
			DeleteGraphPlot(panel, PANEL_GRAPH, -1, VAL_IMMEDIATE_DRAW);
			PlotXY(panel, PANEL_GRAPH, axaX, waveData+start, fs, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
			
			//PlotY(panel, PANEL_GRAPH, waveData+start, fs, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
			SetCtrlVal(panel,PANEL_START, start/fs);
			SetCtrlVal(panel,PANEL_STOP1, (start+fs)/fs);
			
			GetCtrlVal(panel,PANEL_STOP1, &stopc); 
			GetCtrlVal(panel,PANEL_START, &startc);
			
		}
			
			break;
		}
	}
	return 0;
}

int CVICALLBACK filter (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	
	double alpha;
	int WindowDim;
	int optiune;
	
	int *axaX=(int*)calloc(npoints, sizeof(int));
	switch (event)
	{
			case EVENT_COMMIT:
				DeleteGraphPlot(panel, PANEL_GRAPHF, -1, VAL_IMMEDIATE_DRAW);
				GetCtrlVal(panel,PANEL_STOP1, &stopc);
				GetCtrlVal(panel,PANEL_START, &startc);
				GetCtrlVal(panel, PANEL_RING, &optiune);
				if(optiune==1)
				{
					
					double *filt=(double*)calloc(npoints, sizeof(double));
					filt[0]=0;
					GetCtrlVal(panel,PANEL_ALPHA, &alpha);	
					for(int i=1;i<npoints;++i)
					{
						filt[i]=(1-alpha)*filt[i-1]+alpha*waveData[i];
					}
					if(startc==stopc)
					{
						for(int i=0;i<npoints;++i)
							axaX[i]=i;
						PlotXY(panel, PANEL_GRAPHF, axaX,filt, npoints, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);

					}
					else
					{
						for(int i=0;i<fs;++i)
							axaX[i]=i+start;
						PlotXY(panel, PANEL_GRAPHF, axaX, filt+start, fs, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);

					}
					
			
					
				}
				if(optiune==2)
				{
					GetCtrlVal(panel,PANEL_DIM, &WindowDim);	
					
					double *semnalfiltrat=(double*)calloc(npoints, sizeof(double));
					int k=0;
					double s=0;
					for(int i=0;i<npoints-WindowDim;++i)
					{
						if(i==0)
						{
							for(int j=0;j<=WindowDim-1;++j)
								{
									s+=waveData[j];
								}
								semnalfiltrat[k]=s/WindowDim;
								k++;
						}
							
						else		
						{
							s=s-waveData[i-1]+waveData[WindowDim+i-1];
							semnalfiltrat[k]=s/WindowDim;
							k++;
							 
						}
					}
					for(int i=fs-WindowDim;i<fs;++i)
					{
						semnalfiltrat[k]=semnalfiltrat[k-1];  
						k++;
					}
				if(startc==stopc)
					{
						for(int i=0;i<npoints;++i)
							axaX[i]=i;
						PlotXY(panel, PANEL_GRAPHF, axaX,semnalfiltrat, npoints, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
					
					}
					else
					{
						for(int i=0;i<fs;++i)
							axaX[i]=i+start;
						PlotXY(panel, PANEL_GRAPHF, axaX, semnalfiltrat+start, fs, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);

					}	
				}
				
			break;	
	}
		
	return 0;
}


int CVICALLBACK SaveImg (int panel, int control, int event,
						void *callbackData, int eventData1, int eventData2)
{
	int year = 0;
	int month = 0;
	int day = 0;
	int hour = 0;
	int minute = 0;
	int second = 0;
	int imgHandle;  
	char fileName[256] = {0},fileName1[256]={0};
	switch (event)
	{
		case EVENT_COMMIT:
			GetSystemDate(&month, &day, &year);
			GetSystemTime(&hour, &minute, &second);
			
				
			sprintf(fileName, "ROW_DATA_%4d.%2d.%2d_%2d-%2d-%2d.jpg", year, month, day, hour, minute, second);
			
			GetCtrlDisplayBitmap(panel,PANEL_GRAPH,1,&imgHandle);
			SaveBitmapToJPEGFile(imgHandle,fileName,JPEG_PROGRESSIVE,100);
			
			
			sprintf(fileName1, "FILTERED_DATA__%4d.%2d.%2d_%2d-%2d-%2d.jpg", year, month, day, hour, minute, second);
		
			GetCtrlDisplayBitmap(panel,PANEL_GRAPHF,1,&imgHandle);
			SaveBitmapToJPEGFile(imgHandle,fileName1,JPEG_PROGRESSIVE,100);
			
			DiscardBitmap(imgHandle);
			
			break;
	}
	return 0;
}


int CVICALLBACK Derivata (int panel, int control, int event,
							  void *callbackData, int eventData1, int eventData2)
{
	
	int *axaX1=(int*)calloc(npoints, sizeof(int));
	
	
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panel,PANEL_STOP1, &stopc);
			GetCtrlVal(panel,PANEL_START, &startc);
			if (npoints == 0)
				return -1;
			
			deriv = (double *) calloc(npoints, sizeof(double));
			DifferenceEx (waveData, npoints, 1.0, waveData + 1, 1, waveData + npoints - 1, 1, DIFF_SECOND_ORDER_CENTRAL, deriv); 
			if(startc==stopc)
				{
					for(int i=0;i<npoints;++i)
						axaX1[i]=i;
					PlotXY(panel, PANEL_GRAPHF, axaX1, deriv, npoints, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_WHITE);   
				
				}
			else
				{
					for(int i=0;i<fs;++i)
						axaX1[i]=i+start;
					PlotXY(panel, PANEL_GRAPHF, axaX1, deriv+start, fs, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_WHITE);   
					
				}
		
			break;
	}
	return 0;
}


int CVICALLBACK Anvelopa1(int panel, int control, int event,
							  void *callbackData, int eventData1, int eventData2)
{
	
	double* peakLocations;
	
	double* peakAmplitudes;
	double* peakSecondDerivatives;
	double* waveData1;
	double* waveData2;
	
	double *index1;
	double *index2;
	int m;
	ssize_t count;
	
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panel,PANEL_STOP1, &stopc);
			GetCtrlVal(panel,PANEL_START, &startc);
			peakLocations = (double *) calloc(npoints, sizeof(double));
			peakAmplitudes = (double *) calloc(npoints, sizeof(double));
			peakSecondDerivatives = (double *) calloc(npoints, sizeof(double));
			PeakDetector (waveData, npoints, 100, 3, 0, 1, 1, &count, &peakLocations, &peakAmplitudes, &peakSecondDerivatives);	
			double max;
			int maxindex;
			waveData1 = (double *) calloc(npoints, sizeof(double));
			waveData2 = (double *) calloc(npoints, sizeof(double));
			
			index1 = (double *) calloc(npoints, sizeof(double));
			index2 = (double *) calloc(npoints, sizeof(double));
			
			int k=0;
			for(int i=0;i<count;i=i+10)
			{
				maxindex=peakLocations[i];
				max=peakAmplitudes[i];
				for(int j=i;j<i+10;++j)	
				 {
					 if( peakAmplitudes[j]> max)
					 {
							max= peakAmplitudes[j];
							maxindex=peakLocations[j];
					 }
				}waveData1[k]=max;
				index1[k]=maxindex;
				k++;
			}
			
			
			if(startc==stopc)
					{
							PlotXY(panel, PANEL_GRAPH, index1, waveData1, k, VAL_DOUBLE, VAL_DOUBLE, VAL_FAT_LINE, VAL_SOLID_CIRCLE, VAL_SOLID, VAL_SOLID, VAL_WHITE);
							
					}
			else
					{		
							m=0;
							for(int i=0;i<k;++i)
							{
								if(index1[i]>=start && index1[i]<=start+fs)
								{
									
									index2[m]=index1[i];

									waveData2[m]=waveData1[i];
									m++;
								}
									
							}
							PlotXY(panel, PANEL_GRAPH, index2, waveData2, m, VAL_DOUBLE, VAL_DOUBLE, VAL_FAT_LINE, VAL_SOLID_CIRCLE, VAL_SOLID, VAL_SOLID, VAL_WHITE);
						
					}

			break;
}
	return 0;
}




int CVICALLBACK OnTikTak (int panel, int control, int event,
						  void *callbackData, int eventData1, int eventData2)
{
	
	switch (event)
	{
		case EVENT_TIMER_TICK:
			
			break;
	}
	return 0;
}





int CVICALLBACK Filtrare (int panel, int control, int event,
						  void *callbackData, int eventData1, int eventData2)
{
	
	
	switch (event)
	{
		case EVENT_COMMIT:
			DeleteGraphPlot(panel, PANEL_FREQ_GRAPH_5, -1, VAL_IMMEDIATE_DRAW);
			DeleteGraphPlot(panel, PANEL_FREQ_GRAPH_6, -1, VAL_IMMEDIATE_DRAW);
			DeleteGraphPlot(panel, PANEL_FREQ_GRAPH_7, -1, VAL_IMMEDIATE_DRAW);
			DeleteGraphPlot(panel, PANEL_FREQ_GRAPH, -1, VAL_IMMEDIATE_DRAW);
				GetCtrlVal(panelHandleFreq, PANEL_FREQ_RINGSLIDE, &N);
				int* axaX=0;
				int type=0;
				int i, j;
				waveDataD = (double *) calloc(npoints, sizeof(double));
				Copy1D(waveData, npoints ,waveDataD);
   				for (i = 0, j = 0; i < npoints; i +=2, j++) {
        			waveDataD[j] = waveDataD[i];
				}
				fs=fs/2;
				
				FilteredDataa=(double*)calloc(npoints,sizeof(double));
				waveCop=(double*)calloc(N,sizeof(double));
				axaX=(int*)calloc(N,sizeof(int));
				GetCtrlVal(panelHandleFreq, PANEL_FREQ_FilterRing, &type);
				if(type==2)
				{
					GetCtrlVal(panelHandleFreq, PANEL_FREQ_RING, &secunda);
					Copy1D(waveDataD+fs*(secunda-1), N ,waveCop);
					SavitzkyGolayFiltering (waveCop, N, 2, 27, NULL, FilteredDataa);
				}
				
				if(type==1)
				{
					GetCtrlVal(panelHandleFreq, PANEL_FREQ_RING, &secunda);
					int numDen=3;
					int numNum=3;
					Copy1D(waveDataD+fs*(secunda-1), N ,waveCop);
					double* Denominator=(double*)calloc(numDen+1,sizeof(double));
					double* Numerator=(double*)calloc(numNum+1,sizeof(double));
				
					
					DeleteGraphPlot (panel,PANEL_FREQ_GRAPH_6, -1, VAL_IMMEDIATE_DRAW);
				
					Numerator[0]= 0.91941075;//1316481089645549218403175473213
					Numerator[1]=-1.82694368;//5272612416525817025103606283665
					Numerator[2]= 0.91941075;//1316481089645549218403175473213
					Denominator[0]=1; 
					Denominator[1]=-1.82694368;//5272612416525817025103606283665
					Denominator[2]= 0.83882150;//2632962179291098436806350946426
					FilteredDataa[0]=waveCop[0];FilteredDataa[1]=waveCop[1];
					for(int i=2;i<N;i++)
					{
						FilteredDataa[i]=Numerator[0]*waveCop[i]+Numerator[1]*waveCop[i-1]+Numerator[2]*waveCop[i-2]-Denominator[1]*FilteredDataa[i-1]-Denominator[2]*FilteredDataa[i-2];
					}
					
				}
				for(int i=0;i<N;++i)
				{
					axaX[i]=i+fs*(secunda-1);
				}
				DeleteGraphPlot(panel, PANEL_FREQ_GRAPH_4, -1, VAL_IMMEDIATE_DRAW);
				PlotXY(panel, PANEL_FREQ_GRAPH_4, axaX, waveCop, N, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
			
				DeleteGraphPlot(panel, PANEL_FREQ_GRAPH_8, -1, VAL_IMMEDIATE_DRAW);
				PlotXY(panel, PANEL_FREQ_GRAPH_8, axaX, FilteredDataa, N, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
			
			break;
	}
	return 0;
}



int CVICALLBACK Spectru (int panel, int control, int event,
						 void *callbackData, int eventData1, int eventData2)
{	
	int type;
	char unit[32]="V";
	convertedSpectrum=(double*)calloc(N,sizeof(double));
	autoSpectrum=(double*)calloc(N,sizeof(double));
	convertedSpectrumF=(double*)calloc(N,sizeof(double));
	autoSpectrumF=(double*)calloc(N,sizeof(double));
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelHandleFreq, PANEL_FREQ_RING, &secunda);
			GetCtrlVal(panelHandleFreq, PANEL_FREQ_WindowRing, &type);
			if(type==2)
			{
				//spectru nefiltrat
				ScaledWindowEx (waveCop, N, RECTANGLE, 0, &winConst);
				//afisare ferestruire
				DeleteGraphPlot (panel,PANEL_FREQ_GRAPH_6, -1, VAL_IMMEDIATE_DRAW);	
				PlotY (panel,PANEL_FREQ_GRAPH_6, waveCop, N, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
				//
				AutoPowerSpectrum(waveCop,N, 1.0/(sampleRate/2), autoSpectrum, &df);
				PowerFrequencyEstimate(autoSpectrum,N/2,-1.0,winConst,df,7,&freqPeak,&powerPeak);
				
				SetCtrlVal(panelHandleFreq, PANEL_FREQ_FREQPEAK, freqPeak);
				SetCtrlVal(panelHandleFreq, PANEL_FREQ_PEAKPOWER, powerPeak);
				
				SpectrumUnitConversion(autoSpectrum, N/2, 0, SCALING_MODE_LINEAR, DISPLAY_UNIT_VRMS, df, winConst, convertedSpectrum, unit);
				DeleteGraphPlot(panel,PANEL_FREQ_GRAPH,-1,VAL_IMMEDIATE_DRAW);
				PlotWaveform(panel,  PANEL_FREQ_GRAPH, convertedSpectrum, N/2,VAL_DOUBLE, 1.0, 0.0, 0.0, df,VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);	
				
				//spectru filtrat
				ScaledWindowEx (FilteredDataa, N, RECTANGLE, 0, &winConstF);
				//afisare ferestruire
				DeleteGraphPlot (panel,PANEL_FREQ_GRAPH_7, -1, VAL_IMMEDIATE_DRAW);	
				PlotY (panel,PANEL_FREQ_GRAPH_7, FilteredDataa, N, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
				
				AutoPowerSpectrum(FilteredDataa,N, 1.0/(sampleRate/2), autoSpectrumF, &df);
				PowerFrequencyEstimate(autoSpectrumF, N/2,-1.0,winConstF,df,7,&freqPeak,&powerPeak);
				
				//SetCtrlVal(panelHandleFreq, PANEL_FREQ_FREQPEAK, freqPeak);
				//SetCtrlVal(panelHandleFreq, PANEL_FREQ_PEAKPOWER, powerPeak);
				
				SpectrumUnitConversion(autoSpectrumF, N/2, 0, SCALING_MODE_LINEAR, DISPLAY_UNIT_VRMS, df, winConstF, convertedSpectrumF, unit);
				DeleteGraphPlot(panel,PANEL_FREQ_GRAPH_5,-1,VAL_IMMEDIATE_DRAW);
				PlotWaveform(panel,  PANEL_FREQ_GRAPH_5, convertedSpectrumF, N/2,VAL_DOUBLE, 1.0, 0.0, 0.0, df,VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);	
				
				//DeleteGraphPlot (panel,PANEL_FREQ_GRAPH_6, -1, VAL_IMMEDIATE_DRAW);
				//PlotY (panel,PANEL_FREQ_GRAPH_6, winConstF, N, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
						
				}
			
			if(type==1)
			{
				//spectru nefiltrat
				ScaledWindowEx (waveCop, N, GAUSSIAN, 0, &winConst);
				DeleteGraphPlot (panel,PANEL_FREQ_GRAPH_6, -1, VAL_IMMEDIATE_DRAW);	
				PlotY (panel,PANEL_FREQ_GRAPH_6, waveCop, N, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
				
				AutoPowerSpectrum(waveCop,N, 1.0/(sampleRate/2), autoSpectrum, &df);
				PowerFrequencyEstimate(autoSpectrum,N/2,-1.0,winConst,df,7,&freqPeak,&powerPeak);
				
				SetCtrlVal(panelHandleFreq, PANEL_FREQ_FREQPEAK, freqPeak);
				SetCtrlVal(panelHandleFreq, PANEL_FREQ_PEAKPOWER, powerPeak);
				
				SpectrumUnitConversion(autoSpectrum, N/2, 0, SCALING_MODE_LINEAR, DISPLAY_UNIT_VRMS, df, winConst, convertedSpectrum, unit);
				DeleteGraphPlot(panel,PANEL_FREQ_GRAPH,-1,VAL_IMMEDIATE_DRAW);
				PlotWaveform(panel,  PANEL_FREQ_GRAPH, convertedSpectrum, N/2,VAL_DOUBLE, 1.0, 0.0, 0.0, df,VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);	
				
				//spectru filtrat
				ScaledWindowEx (FilteredDataa, N, GAUSSIAN, 0, &winConstF);
				
				DeleteGraphPlot (panel,PANEL_FREQ_GRAPH_7, -1, VAL_IMMEDIATE_DRAW);	
				PlotY (panel,PANEL_FREQ_GRAPH_7, FilteredDataa, N, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);
				
				AutoPowerSpectrum(FilteredDataa,N, 1.0/(sampleRate/2), autoSpectrumF, &df);
				PowerFrequencyEstimate(autoSpectrumF, N/2,-1.0,winConstF,df,7,&freqPeak,&powerPeak);
				
				//SetCtrlVal(panelHandleFreq, PANEL_FREQ_FREQPEAK, freqPeak);
				//SetCtrlVal(panelHandleFreq, PANEL_FREQ_PEAKPOWER, powerPeak);
				
				SpectrumUnitConversion(autoSpectrumF, N/2, 0, SCALING_MODE_LINEAR, DISPLAY_UNIT_VRMS, df, winConstF, convertedSpectrumF, unit);
				DeleteGraphPlot(panel,PANEL_FREQ_GRAPH_5,-1,VAL_IMMEDIATE_DRAW);
				PlotWaveform(panel,  PANEL_FREQ_GRAPH_5, convertedSpectrumF, N/2,VAL_DOUBLE, 1.0, 0.0, 0.0, df,VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_CYAN);	
			}
				
			
			break;
	}
	return 0;

}


int CVICALLBACK SaveFrecv (int panel, int control, int event,
						   void *callbackData, int eventData1, int eventData2)
{
	int year = 0;
	int month = 0;
	int day = 0;
	int hour = 0;
	int minute = 0;
	int second = 0;
	int imgHandle;  
	char fileName[256] = {0},fileName1[256]={0};
	switch (event)
	{
		case EVENT_COMMIT:
			GetSystemDate(&month, &day, &year);
			GetSystemTime(&hour, &minute, &second);
			
				
			sprintf(fileName, "ROW_DATA_Frecv__%4d.%2d.%2d_%2d-%2d-%2d.jpg", year, month, day, hour, minute, second);
			
			GetCtrlDisplayBitmap(panel,PANEL_FREQ_GRAPH_4,1,&imgHandle);
			SaveBitmapToJPEGFile(imgHandle,fileName,JPEG_PROGRESSIVE,100);
			
			
			sprintf(fileName1, "FILTERED_DATA_Frecv__%4d.%2d.%2d_%2d-%2d-%2d.jpg", year, month, day, hour, minute, second);
		
			GetCtrlDisplayBitmap(panel,PANEL_FREQ_GRAPH_8,1,&imgHandle);
			SaveBitmapToJPEGFile(imgHandle,fileName1,JPEG_PROGRESSIVE,100);
			
			
			sprintf(fileName1, "WINDOW_Frecv__%4d.%2d.%2d_%2d-%2d-%2d.jpg", year, month, day, hour, minute, second);
		
			GetCtrlDisplayBitmap(panel,PANEL_FREQ_GRAPH_6,1,&imgHandle);
			SaveBitmapToJPEGFile(imgHandle,fileName1,JPEG_PROGRESSIVE,100);
			
			
			sprintf(fileName1, "FILTERED_WINDOW_Frecv__%4d.%2d.%2d_%2d-%2d-%2d.jpg", year, month, day, hour, minute, second);
		
			GetCtrlDisplayBitmap(panel,PANEL_FREQ_GRAPH_7,1,&imgHandle);
			SaveBitmapToJPEGFile(imgHandle,fileName1,JPEG_PROGRESSIVE,100);
			
			
			sprintf(fileName1, "SPECTRUM__%4d.%2d.%2d_%2d-%2d-%2d.jpg", year, month, day, hour, minute, second);
		
		 	GetCtrlDisplayBitmap(panel,PANEL_FREQ_GRAPH_5,1,&imgHandle);
			SaveBitmapToJPEGFile(imgHandle,fileName1,JPEG_PROGRESSIVE,100);
			
			
			sprintf(fileName1, "FILTERED_SPECTRUM__%4d.%2d.%2d_%2d-%2d-%2d.jpg", year, month, day, hour, minute, second);
		
			GetCtrlDisplayBitmap(panel,PANEL_FREQ_GRAPH,1,&imgHandle);
			SaveBitmapToJPEGFile(imgHandle,fileName1,JPEG_PROGRESSIVE,100);
			DiscardBitmap(imgHandle);
	
			break;
	}
	return 0;
}



