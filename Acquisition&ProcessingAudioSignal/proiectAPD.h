/**************************************************************************/
/* LabWindows/CVI User Interface Resource (UIR) Include File              */
/*                                                                        */
/* WARNING: Do not add to, delete from, or otherwise modify the contents  */
/*          of this include file.                                         */
/**************************************************************************/

#include <userint.h>

#ifdef __cplusplus
    extern "C" {
#endif

     /* Panels and Controls: */

#define  PANEL                            1       /* callback function: On_Panel */
#define  PANEL_GRAPHF                     2       /* control type: graph, callback function: (none) */
#define  PANEL_GRAPH                      3       /* control type: graph, callback function: (none) */
#define  PANEL_COMMANDBUTTON              4       /* control type: command, callback function: On_Prev */
#define  PANEL_NEXT                       5       /* control type: command, callback function: On_Next */
#define  PANEL_COMMANDBUTTON_3            6       /* control type: command, callback function: load */
#define  PANEL_RING                       7       /* control type: ring, callback function: (none) */
#define  PANEL_FIILTRU                    8       /* control type: command, callback function: filter */
#define  PANEL_MININDEX                   9       /* control type: numeric, callback function: (none) */
#define  PANEL_MIN                        10      /* control type: numeric, callback function: (none) */
#define  PANEL_MAXINDEX                   11      /* control type: numeric, callback function: (none) */
#define  PANEL_MAX                        12      /* control type: numeric, callback function: (none) */
#define  PANEL_DISP                       13      /* control type: numeric, callback function: (none) */
#define  PANEL_MEDIANA                    14      /* control type: numeric, callback function: (none) */
#define  PANEL_MEDIE                      15      /* control type: numeric, callback function: (none) */
#define  PANEL_GRAPH_3                    16      /* control type: graph, callback function: (none) */
#define  PANEL_ZERO                       17      /* control type: numeric, callback function: (none) */
#define  PANEL_STOP1                      18      /* control type: numeric, callback function: (none) */
#define  PANEL_START                      19      /* control type: numeric, callback function: (none) */
#define  PANEL_ALPHA                      20      /* control type: numeric, callback function: (none) */
#define  PANEL_DIM                        21      /* control type: numeric, callback function: (none) */
#define  PANEL_KTS                        22      /* control type: numeric, callback function: (none) */
#define  PANEL_SK                         23      /* control type: numeric, callback function: (none) */
#define  PANEL_COMMANDBUTTON_2            24      /* control type: command, callback function: SaveImg */
#define  PANEL_DERIV                      25      /* control type: command, callback function: Derivata */
#define  PANEL_COMMANDBUTTON_4            26      /* control type: command, callback function: Anvelopa1 */
#define  PANEL_SPLITTER_20                27      /* control type: splitter, callback function: (none) */
#define  PANEL_SPLITTER_16                28      /* control type: splitter, callback function: (none) */
#define  PANEL_SPLITTER_2                 29      /* control type: splitter, callback function: (none) */
#define  PANEL_SPLITTER                   30      /* control type: splitter, callback function: (none) */
#define  PANEL_SPLITTER_17                31      /* control type: splitter, callback function: (none) */
#define  PANEL_SPLITTER_19                32      /* control type: splitter, callback function: (none) */
#define  PANEL_SPLITTER_21                33      /* control type: splitter, callback function: (none) */
#define  PANEL_SPLITTER_18                34      /* control type: splitter, callback function: (none) */
#define  PANEL_SPLITTER_7                 35      /* control type: splitter, callback function: (none) */
#define  PANEL_BINARYSWITCH               36      /* control type: binary, callback function: OnSwitchPanelCB */

#define  PANEL_FREQ                       2       /* callback function: On_Panel_Freq */
#define  PANEL_FREQ_BINARYSWITCH          2       /* control type: binary, callback function: OnSwitchPanelCB */
#define  PANEL_FREQ_GRAPH_5               3       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_GRAPH_7               4       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_GRAPH_6               5       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_GRAPH_8               6       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_GRAPH_4               7       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_GRAPH                 8       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_RINGSLIDE             9       /* control type: slide, callback function: (none) */
#define  PANEL_FREQ_SPECTRU               10      /* control type: command, callback function: Spectru */
#define  PANEL_FREQ_PEAKPOWER             11      /* control type: numeric, callback function: (none) */
#define  PANEL_FREQ_FREQPEAK              12      /* control type: numeric, callback function: (none) */
#define  PANEL_FREQ_RING                  13      /* control type: ring, callback function: (none) */
#define  PANEL_FREQ_TIMER                 14      /* control type: timer, callback function: OnTikTak */
#define  PANEL_FREQ_WindowRing            15      /* control type: ring, callback function: (none) */
#define  PANEL_FREQ_FilterRing            16      /* control type: ring, callback function: (none) */
#define  PANEL_FREQ_Filtrare              17      /* control type: command, callback function: Filtrare */
#define  PANEL_FREQ_save                  18      /* control type: command, callback function: SaveFrecv */


     /* Control Arrays: */

          /* (no control arrays in the resource file) */


     /* Menu Bars, Menus, and Menu Items: */

#define  MENUBAR                          1
#define  MENUBAR_Proiect_APD              2


     /* Callback Prototypes: */

int  CVICALLBACK Anvelopa1(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Derivata(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK filter(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Filtrare(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK load(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK On_Next(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK On_Panel(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK On_Panel_Freq(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK On_Prev(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnSwitchPanelCB(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnTikTak(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SaveFrecv(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SaveImg(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Spectru(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif