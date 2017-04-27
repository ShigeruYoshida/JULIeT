package iceCube.uhe.analysis;

import hep.aida.*;

public class JaidaPlotStyleSetter {

    public static void setPlotStyle(IPlotterStyle style,
				       String xLabel, String yLabel){

	style.xAxisStyle().setLabel(xLabel);
	style.yAxisStyle().setLabel(yLabel);
	style.xAxisStyle().tickLabelStyle().setBold(true);
	style.xAxisStyle().tickLabelStyle().setItalic(true);
	style.xAxisStyle().tickLabelStyle().setFontSize(14);
	style.yAxisStyle().tickLabelStyle().setBold(true);
	style.yAxisStyle().tickLabelStyle().setFontSize(14);
	style.titleStyle().textStyle().setFontSize(20);
    }

    public static void setPlotStyle(IPlotterStyle style,
				    String xLabel, String yLabel,
				    String parameterName, String option){

	setPlotStyle(style,xLabel,yLabel);
	style.setParameter(parameterName,option);
    }
}

