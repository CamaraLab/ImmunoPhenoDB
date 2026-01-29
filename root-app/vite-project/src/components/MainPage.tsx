import React, { useCallback, useMemo, useState, useEffect, useRef } from 'react'
import TissuesTable from './TissuesTable'
import RadioSearchField from './RadioSearchField';
import SelectionTable from './SelectionTable';
import Cart from './Cart';
import ModelResultsTable from './ModelResultsTable';
import ModelResultsPlot from './ModelResultsPlot';

import {
  AppBar,
  Toolbar,
  Typography,
  Box,
  CssBaseline,
  Grid, // Import Grid
  Stack, // Stack is useful for vertical arrangement within columns
  Paper, // Use Paper for simple text sections too
  Link
} from '@mui/material';

const MainPage = () => {

  // Props for TissueTable.tsx
  const [tissueRows, setTissueRows] = useState({});

  // Note:
  // When radioValue = antibodies, we are calling find_celltypes (given antibodies)
  // When radioValue = cell types, we are calling find_antibodies (given cell types)

  // Props for RadioSearchField.tsx
  const [radioValue, setRadioValue] = useState("Antibodies");   // Radio button selection
  const [searchQuery, setSearchQuery] = useState("");           // Textfield query
  const [searchResults, setSearchResults] = useState([]);       // Rows in the selection table
  const [searchColumns, setSearchColumns] = useState<{ accessorKey: string; header: string }[]>([]);       // Column names in selection table
  const [searchIsLoading, setSearchIsLoading] = useState(false);// Loading animation for selection table
  const [cartColumns, setCartColumns] = useState<{ accessorKey: string; header: string }[]>([]);            // Set columns for cart later


  // Not sure what these do atm -- now i see, it's for later. when you change the radio, clear these values too
  // These are used later for the modelResult table, holding the normalized values
  const [modelResults, setModelResults] = useState([]);         // Model results table (used for cell types radioValue option)
  const [modelColumns, setModelColumns] = useState([]);         // Column names for result table
  const [modelIsLoading, setModelIsLoading] = useState(false);  // Loading animation for result table
  
  const [modelAbKeys, setModelAbKeys] = useState([])   // what? -- now i see, this api returns a dictionary where keys are ab names. Linear Mixed Model Antibody Keys (for antibodies radiovalue/find_celltypes)
  const [modelResultsDict, setModelResultsDict] = useState({})
  console.log("modelResults:", modelResults)

  // Props for SelectionTable.tsx
  const [searchRowSelection, setSearchRowSelection] = useState([]) // Starting off with nothing selected
  const [cartData, setCartData] = useState([])
  const [emptyResultsError, setEmptyResultsError] = useState(false) // Error handling for table AND plots

  // Props for Cart.tsx use the ones defined above
  const [modelPlotIsLoading, setModelPlotIsLoading] = useState(false) // load plot animation
  const [modelPlotResults, setModelPlotResults] = useState([]) // needs to have empty array?, or else we read undefined ((used for cell types radioValue option))
  console.log("setModelPlotResults", modelPlotResults)
//   const [modelPlotResultsDict, setModelPlotResultsDict] = useState([]) // used for antibodies radioValue
  const [modelPlotResultsDict, setModelPlotResultsDict] = useState(new Map()); // using Map



  // Props for ModelResultsTable.tsx
  const [tableIndex, setTableIndex] = useState(0) // Used for "Antibodies"

  // Interfaces for format of selection table, and format of search cart
interface AntibodyTable {
    idAntibody: string;
    abName: string;
    abTarget: string;
    clonality: string;
    citation: string;
    comments: string;
    cloneID: string;
    host: string;
    vendor: string;
    catalogNum: string;
    idExperiment: string[];
  };

interface CellTypeTable {
    idCL: string;
    label: string;
    idExperiment: string[];
}

interface AntibodyCart {
    idAntibody: string;
    abTarget: string;
}

interface CellTypeCart {
    idCL: string;
    label: string;
}

  const AntibodyColumns = [
      {
        accessorKey: 'idAntibody',
        header: 'Antibody ID',
        // size: 125,
      },
      {
        accessorKey: 'abName',
        header: 'Antibody Name',
      },
      {
        accessorKey: 'abTarget',
        header: 'Antibody Target',
      },
      {
        accessorKey: 'clonality',
        header: 'Clonality',
      },
      {
        accessorKey: 'citation',
        header: 'Citation',
      },
      {
        accessorKey: 'comments',
        header: 'Comments',
      },
      {
        accessorKey: 'cloneID',
        header: 'Clone ID',
      },
      {
        accessorKey: 'host',
        header: 'Host Organism',
      },
      {
        accessorKey: 'vendor',
        header: 'Vendor',
      },
      {
        accessorKey: 'catalogNum',
        header: 'Catalog Number',
      },
      {
        accessorKey: 'idExperiment',
        header: 'Found in Experiment(s)',
      },
      //end
    ]

  const CellTypeColumns = [
        {
          accessorKey: "idCL",
          header: "Cell Ontology ID",
        },
        {
          accessorKey: "label",
          header: "Cell Type"
        },
        {
          accessorKey: 'idExperiment',
          header: 'Found in Experiment(s)',
        }
        //end
    ]

  const AntibodyCart = [
        {
            accessorKey: 'idAntibody',
            header: 'Antibody ID',
            // size: 125,
        },
        {
            accessorKey: 'abTarget',
            header: 'Antibody Target',
        },
      ]

  const CellTypeCart = [
        {
            accessorKey: "idCL",
            header: "Cell Ontology ID",
        },
        {
            accessorKey: "label",
            header: "Cell Type"
        },
      ]

// // // Need to add useEffect hook here to reset table, plot, cart values when radio changes?
// // // On initial load, display antibodies table and cart
// useEffect(() => {
//     setSearchColumns(AntibodyColumns)
//     setCartColumns(AntibodyCart)
// }, [])

// // Clear any text and search table results upon load/reload
// useEffect(() => {
//     setSearchQuery("")
//     setSearchResults([])
//     setSearchIsLoading(false)
//     if (radioValue === "Antibodies") {
//         setSearchColumns(AntibodyColumns)
//         setCartColumns(AntibodyCart)
        
//       } else if (radioValue === "Cell Types") {
//         setSearchColumns(CellTypeColumns)
//         setCartColumns(CellTypeCart)
//       }
//  }, [radioValue]) //dependency array, useEffect will execute whatever when this array changes
const appBarHeight: number = 64; // Standard AppBar height

  return (
    <div>
        {/* <div>
            <TissuesTable 
              setTissueRows={setTissueRows}
            />
        </div>
        <br></br>
        <div>
            <RadioSearchField 
              radioValue={radioValue}
              setRadioValue={setRadioValue}
              setSearchQuery={setSearchQuery}
              setSearchResults={setSearchResults}
              setSearchColumns={setSearchColumns}
              setSearchIsLoading={setSearchIsLoading}
              setCartColumns={setCartColumns}
              setModelColumns={setModelColumns}
              setModelResults={setModelResults}
              setModelAbKeys={setModelAbKeys}

              setModelResultsDict={setModelResultsDict}
              setModelPlotResults={setModelPlotResults}
              setModelPlotResultsDict={setModelPlotResultsDict}
            />
        </div>
        <br></br>
        <div>
            <SelectionTable 
              searchResults={searchResults}
              searchColumns={searchColumns}
              searchIsLoading={searchIsLoading}
              searchRowSelection={searchRowSelection}
              setSearchRowSelection={setSearchRowSelection}
              setCartData={setCartData}
              radioValue={radioValue}
            />
        </div>
        <div>
            <Cart 
              cartData={cartData}
              cartColumns={cartColumns}
              tissueRows={tissueRows}
              setCartData={setCartData}
              setSearchRowSelection={setSearchRowSelection}
              radioValue={radioValue}

              setEmptyResultsError={setEmptyResultsError}

              setModelColumns={setModelColumns}
              setModelIsLoading={setModelIsLoading}
              setModelPlotIsLoading={setModelPlotIsLoading}
              
              setModelResults={setModelResults}
              setModelPlotResults={setModelPlotResults}

              setModelResultsDict={setModelResultsDict}
              setModelPlotResultsDict={setModelPlotResultsDict}
              setModelAbKeys={setModelAbKeys}
              
            />
        </div>
        <br></br>
        <div>
            <ModelResultsTable 
              radioValue={radioValue}
              modelColumns={modelColumns}
              modelResults={modelResults}
              modelIsLoading={modelIsLoading}
              
              emptyResultsError={emptyResultsError}
              setEmptyResultsError={setEmptyResultsError}

              modelAbKeys={modelAbKeys}
              tableIndex={tableIndex}
              setTableIndex={setTableIndex}
              modelResultsDict={modelResultsDict}
            />
        </div>
        <div>
            <ModelResultsPlot
              radioValue={radioValue} 
              modelPlotResults={modelPlotResults}
              modelPlotIsLoading={modelPlotIsLoading}

              modelAbKeys={modelAbKeys}
              tableIndex={tableIndex}
              modelPlotResultsDict={modelPlotResultsDict}
              emptyResultsError={emptyResultsError}
            />
        </div>
        <br></br> */}
        

        {/* Experimental Below */}
        <Box sx={{ display: 'flex', flexDirection: 'column', minHeight: '100vh' }}>
        <CssBaseline />

      {/* --- Navigation Bar (NavBar.tsx equivalent) IGNORE for now --- */}
      

      {/* --- Main Content Area (Below AppBar) --- */}
      <Box
        component="main"
        sx={{
          flexGrow: 2, // Takes remaining vertical space
          pt: `${appBarHeight}px`, // Start content below AppBar (padding-top)
          p: 3 // Add overall padding around the grid content
        }}
      >
        {/* --- The Grid Container --- */}
        {/* 'container' defines this Box as the grid container */}
        {/* 'spacing' adds space *between* grid items */}
        <Grid container spacing={3} sx={{justifyContent: "center"}}>
          {/* --- Left Column (Querying) --- */}
          {/* 'item' defines this as a grid item */}
          {/* 'xs={12}' makes it full width on extra-small screens (stacking) */}
          {/* 'md={6}' makes it half width on medium screens and up */}
          {/* Ignore above and use a static width for the left column Grid, which is a Grid item */}
          {/* <Grid item  style={{width: "700px"}}> */}
          <Grid item >
             {/* Use Stack for vertical arrangement within this column */}
            <Stack spacing={3}> {/* Add vertical spacing between elements in this column */}
              {/* <Typography variant="h5" component="h2" gutterBottom>
                Step 1: Specify Criteria
              </Typography> */}
              {/* <TissuesTablePlaceholder /> */}
              <TissuesTable 
              setTissueRows={setTissueRows}
              />
              {/* <CartPlaceholder /> */}
              <Cart 
              cartData={cartData}
              cartColumns={cartColumns}
              tissueRows={tissueRows}
              setCartData={setCartData}
              setSearchRowSelection={setSearchRowSelection}
              radioValue={radioValue}

              setEmptyResultsError={setEmptyResultsError}

              setModelColumns={setModelColumns}
              setModelIsLoading={setModelIsLoading}
              setModelPlotIsLoading={setModelPlotIsLoading}
              
              setModelResults={setModelResults}
              setModelPlotResults={setModelPlotResults}

              setModelResultsDict={setModelResultsDict}
              setModelPlotResultsDict={setModelPlotResultsDict}
              setModelAbKeys={setModelAbKeys}
              
            />
            </Stack>
          </Grid>

          {/* --- Right Column (Results) --- */}
          <Grid item xs={14} md={6}>
            <Stack spacing={3}> {/* Add vertical spacing between elements in this column */}
               {/* <Typography variant="h5" component="h2" gutterBottom>
                 Title: Results
               </Typography> */}
               {/* <RadioSearchFieldPlaceholder /> */}

               <Paper elevation={1} sx={{ p: 2, backgroundColor: '#f5f5f5' }}>
                  <Typography variant="body1">
                  ImmunoPheno is a computational platform for the design and analysis of immunophenotyping experiments based on existing proteotranscriptomic datasets (CITE-seq, Abseq, REAP-seq). 
                  The ImmunoPheno database is a relational system containing structured and harmonized data from published reference proteotranscriptomic studies.
                  This web application allows the user to perform simple queries without needing to install or use the <Link href="https://github.com/CamaraLab/ImmunoPheno">ImmunoPheno</Link> Python package. 
                  It enables the identification of antibodies tailored to marking specific cell populations, assisting in the design of immunophenotyping experiments.
                  To begin querying the databse, select either "Antibodies" or "Cell Types" and a search phrase.
                    </Typography>
               </Paper>
               <RadioSearchField 
              radioValue={radioValue}
              setRadioValue={setRadioValue}
              setSearchQuery={setSearchQuery}
              setSearchResults={setSearchResults}
              setSearchColumns={setSearchColumns}
              setSearchIsLoading={setSearchIsLoading}
              setCartColumns={setCartColumns}
              setModelColumns={setModelColumns}
              setModelResults={setModelResults}
              setModelAbKeys={setModelAbKeys}

              setModelResultsDict={setModelResultsDict}
              setModelPlotResults={setModelPlotResults}
              setModelPlotResultsDict={setModelPlotResultsDict}
            />

               {/* <SelectionTablePlaceholder /> */}
               {/* Text Explanation 1 */}
               {/* <ExplanationPlaceholder text="Text: Explain what the values in the Result table mean" /> */}
               {/* Or use Paper directly: */}
               <Paper elevation={1} sx={{ p: 2, backgroundColor: '#f5f5f5' }}>
                  <Typography variant="body1">
                  This web application implements a linear mixed effects model to identify antibodies that mark one or more cell populations (when searching by cell type) 
                  or to identify cell types marked by one or more antibodies (when searching by antibody).
                  To begin, click “Add” to move selected rows to the menu on the left. When setting up your query, at least one tissue must be selected in the left-hand menu. 
                  Once your selections are made, click “Run” to execute the linear mixed effects model.
                    </Typography>
               </Paper>
               <SelectionTable 
              searchResults={searchResults}
              searchColumns={searchColumns}
              searchIsLoading={searchIsLoading}
              searchRowSelection={searchRowSelection}
              setSearchRowSelection={setSearchRowSelection}
              setCartData={setCartData}
              radioValue={radioValue}
            />

               {/* Text Explanation 1 */}
               {/* <ExplanationPlaceholder text="Text: Explain what the values in the Result table mean" /> */}
               {/* Or use Paper directly: */}
               <Paper elevation={1} sx={{ p: 2, backgroundColor: '#f5f5f5' }}>
                  <Typography variant="body1">
                    When searching for cell types, antibodies with positive coefficients represent those that are upregulated 
                    in the cell populations selected in the menu, while those with negative coefficients are downregulated. 
                    Likewise, when searching for antibodies, cell types with positive coefficients represent those that are marked by the selected antibodies.
                    </Typography>
               </Paper>

               {/* <ModelResultsTablePlaceholder /> */}
               <ModelResultsTable 
              radioValue={radioValue}
              modelColumns={modelColumns}
              modelResults={modelResults}
              modelIsLoading={modelIsLoading}
              
              emptyResultsError={emptyResultsError}
              setEmptyResultsError={setEmptyResultsError}

              modelAbKeys={modelAbKeys}
              tableIndex={tableIndex}
              setTableIndex={setTableIndex}
              modelResultsDict={modelResultsDict}
            />


               {/* Text Explanation 2 */}
               {/* <ExplanationPlaceholder text="Text: Explain what the plots in the Result graph mean" /> */}
               <Paper elevation={1} sx={{ p: 2, backgroundColor: '#f5f5f5' }}>
                  <Typography variant="body1">
                  The distribution of normalized antibody expression values is shown in the plot below. 
                  Antibodies with greater expression in certain cell populations are expected to have higher normalized values. 
                  The plot is sorted in descending order, from highest to lowest expression levels.
                  </Typography>
               </Paper>

               {/* <ModelResultsPlotPlaceholder /> */}
               <ModelResultsPlot
              radioValue={radioValue} 
              modelPlotResults={modelPlotResults}
              modelPlotIsLoading={modelPlotIsLoading}

              modelAbKeys={modelAbKeys}
              tableIndex={tableIndex}
              modelPlotResultsDict={modelPlotResultsDict}
              emptyResultsError={emptyResultsError}
            />
            </Stack>
          </Grid>

        </Grid> {/* End of Grid Container */}
      </Box> {/* End of Main Content Area */}
    </Box>

    </div>
  )
}

export default MainPage