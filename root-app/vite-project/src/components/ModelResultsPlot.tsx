import React, { useCallback, useState, useEffect, memo, useMemo } from 'react'

  
import _ from 'underscore';
import Plotly from 'plotly.js-dist';
import CSS from 'csstype'

// import './css/plotly.css'
import PlotlyLoad from './PlotlyLoad';

interface SummaryStatistics_antibodies {
    idAntibody: string;
    mean: number;
    q1: number;
    median: number;
    q3: number;
    min: number;
    max: number;
  }
  

interface BoxplotProps_antibodies {
    data: SummaryStatistics_antibodies[];
}

const Boxplot_antibodies = memo(({ data }: BoxplotProps_antibodies) => {
    useEffect(() => {
      const plotData = data.map((item) => ({
        type: 'box',
        name: item.idAntibody,
        x: [item.idAntibody], // Position on the x-axis
        q1: [item.q1],
        median: [item.median],
        q3: [item.q3],
        lowerfence: [item.min],
        upperfence: [item.max],
        mean: item.mean,
        boxpoints: false, // Hide individual points
        hoverinfo: 'all',
      }));
  
      const layout = {
        title: 'Antibodies Boxplot',
        xaxis: { 
          title: 'Antibody',
          type: 'category',
          automargin: true, 
        },
        yaxis: { 
          title: 'Normalized Values',
          automargin: true, 
        },
        width: 1300,
        height: 700,
      };
  
      Plotly.newPlot('boxplot_antibodies', plotData as any, layout);
    }, [data]);
  
    return <div id="boxplot_antibodies" />;
  });

interface SummaryStatistics_celltypes {
    idCL: string;
    mean: number;
    q1: number;
    median: number;
    q3: number;
    min: number;
    max: number;
  }
  

interface BoxplotProps_celltypes {
    data: SummaryStatistics_celltypes[];
}

const Boxplot_celltypes = memo(({ data }: BoxplotProps_celltypes) => {
    useEffect(() => {
      const plotData = data.map((item) => ({
        type: 'box',
        name: item.idCL,
        x: [item.idCL], // Position on the x-axis
        q1: [item.q1],
        median: [item.median],
        q3: [item.q3],
        lowerfence: [item.min],
        upperfence: [item.max],
        mean: item.mean,
        boxpoints: false, // Hide individual points
        hoverinfo: 'all',
      }));
  
      const layout = {
        title: 'Cell Types Boxplot',
        xaxis: { 
          title: 'Cell Type',
          type: 'category',
          automargin: true, 
        },
        yaxis: { 
          title: 'Normalized Values',
          automargin: true, 
        },
        width: 1500,
        height: 700,
      };
  
      Plotly.newPlot('boxplot_celltypes', plotData as any, layout);
    }, [data]);
  
    return <div id="boxplot_celltypes" />;
  });



const ModelResultsPlot = ({radioValue,
                           modelPlotResults,
                           modelPlotIsLoading,
                           
                           modelPlotResultsDict,
                           modelAbKeys,
                           tableIndex,

                           emptyResultsError
                          }:
                          {radioValue: any,
                           modelPlotResults: any,
                           modelPlotIsLoading: any,
                           modelPlotResultsDict: any,
                           modelAbKeys: any,
                           tableIndex: any,

                           emptyResultsError: any
                          }) => {
//   let plotData: any = []
  
//   if (modelAbKeys.length > 0) {
//     const currentAntibody = modelAbKeys[tableIndex]
//     console.log("We're using antibody data. Current antibody:", currentAntibody)
//     plotData = modelPlotResultsDict[currentAntibody]
//     if (plotData == null) {
//         plotData = [] // If there were no results to begin with, but we kept modelAbKeys to display an error - set results to empty array
//     }
//     }
//   else {
//     console.log("We're using cell type data")
//     plotData = modelPlotResults
//   }
//   console.log("plotData for plot:", plotData)

const plotData = useMemo(() => {
    if (radioValue === "Antibodies" && !modelPlotIsLoading && modelPlotResultsDict.size > 0) {
      const currentAntibody = modelAbKeys[tableIndex];
    //   return modelPlotResultsDict[currentAntibody] || [];
      return modelPlotResultsDict.get(currentAntibody)
    } else if (radioValue === "Cell Types" && !modelPlotIsLoading && modelPlotResults.length > 0){
      return modelPlotResults
    } else {
      return [];
    }
  }, [modelAbKeys, tableIndex, modelPlotResultsDict, modelPlotResults]);


console.log("modelPlotResults:", modelPlotResults)
console.log("modelPlotResults.length:", modelPlotResults.length)
console.log("what's plotData:", plotData)

const plotData_comp = useMemo(() => {
    if (radioValue === "Antibodies"  && !modelPlotIsLoading && modelPlotResultsDict.size > 0) {
      console.log("triggered plot1?")

      return <Boxplot_celltypes data={plotData} />
    } else if (radioValue === "Cell Types" && !modelPlotIsLoading && modelPlotResults.length > 0) {
      console.log("triggered plot2? HOW")

      return <Boxplot_antibodies data={plotData} />
    } else if (modelPlotIsLoading) {
      return <PlotlyLoad />
    } else {
      return <b></b>
    }
  }, [modelAbKeys, tableIndex, modelPlotResultsDict, modelPlotResults, modelPlotIsLoading]);

  console.log(modelPlotIsLoading)
      
  return (
    <div>
      {/* <div>
        <Boxplot_antibodies data={plotData} />
      </div>
        
      <div>
       <Boxplot_celltypes data={plotData}/>
      </div> */}
      {plotData_comp}
        {/* {!modelPlotIsLoading && modelPlotResultsDict.size == 0 ? <Boxplot_celltypes data={plotData}/> : undefined} */}
    </div>
  )
}

export default ModelResultsPlot