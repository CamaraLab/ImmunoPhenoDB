import React, { useMemo, useState, useEffect, useCallback } from 'react'
import Radio from '@mui/material/Radio';
import RadioGroup from '@mui/material/RadioGroup';
import FormControlLabel from '@mui/material/FormControlLabel';
import FormControl from '@mui/material/FormControl';
import FormLabel from '@mui/material/FormLabel'
// import Card from '@mui/material/Card';
// import Card from 'react-bootstrap/Card'  
import Card from '@mui/material/Card';
import CardContent from '@mui/material/CardContent';
import CardHeader from '@mui/material/CardHeader';
import Typography from '@mui/material/Typography';

import { default as axios } from "axios";

import {
    type MRT_ColumnDef,
} from 'material-react-table';

import { 
    Button, 
    TextField,
} from '@mui/material';


/* 
The RadioSearchField.tsx component consists of two parts: a radio button to indicate 
"Cell Types" or "Antibodies" searching, and a search bar where you type in a query. 
A "Search" button is present to perform the API call.

The results of the API call should follow the format of AntibodyColumns or CellTypeColumns 
for the selection table. Likewise, the cart should have the corresponding AntibodyCart or 
CellTypeCart columns.


*/

// const URL = "http://127.0.0.1:8080"
// const URL = "http://host.docker.internal" // local development

// Deploying on the cloud
// const URL = "http://node:8080" doesn't work
// const URL = "http://ec2-100-28-86-22.compute-1.amazonaws.com" // causes CORS error
// const URL = ""

// Cloud deployment -- new approach (04/17/2025) using react bundled in NGINX
const URL = "" // it should be empty so NGINX intercepts it


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

// Helper functions 
function convertToArrayObject(dictionary: any) {
    const resultspt1: any = []
    for (const [key, values] of Object.entries(dictionary)) {
        for (const [index, value] of Object.entries(values as any[])) {
            (resultspt1[index]  ??= {})[key] = value
        }
    }
    return resultspt1
} 

function includeID(transposeCol: any, keyName:string) {
    const resultspt2: any = []
    for (const [key, values] of Object.entries(transposeCol)) {
        let temp: any = values
        temp[keyName] = key
        resultspt2.push(temp)
    }
    return resultspt2
}

// RadioButton
const RadioSearchField = ({radioValue,
                        setRadioValue, 
                        setSearchQuery, 
                        setSearchResults,
                        setSearchColumns,
                        setSearchIsLoading,
                        setCartColumns,
                        setModelColumns,
                        setModelResults,
                        setModelAbKeys,
                        setModelResultsDict,
                        setModelPlotResults,
                        setModelPlotResultsDict
                        }: 
                        {radioValue: any,
                        setRadioValue: any, 
                        setSearchQuery: any, 
                        setSearchResults: any,
                        setSearchColumns: any,
                        setSearchIsLoading: any,
                        setCartColumns: any,
                        setModelColumns: any,
                        setModelResults: any,
                        setModelAbKeys: any,
                        setModelResultsDict: any,
                        setModelPlotResults: any,
                        setModelPlotResultsDict: any}) => {
    // const [selectedValue, setSelectedValue] = useState("Antibodies"); // No local for this
    const [textInput, setTextInput] = useState(""); // Local, setSearchQuery is for passing props
    const [errorField, setErrorField] = useState(false); // For empty text box
    const [runDisabled, setRunDisabled] = useState(false)

    // useMemo is a hook, needs to be called inside the component
    // Column layout for results table
    const AntibodyColumns = useMemo<MRT_ColumnDef<AntibodyTable>[]>(
      () => [
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
      ],
      [],
    );

    const CellTypeColumns = useMemo<MRT_ColumnDef<CellTypeTable>[]>(
      () => [
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
          },
          //end
      ],
      [],
    );

    // Column layout for cart
    const AntibodyCart = useMemo<MRT_ColumnDef<AntibodyCart>[]>(
      () => [
          {
              accessorKey: 'idAntibody',
              header: 'Antibody ID',
              size: 30,
              maxSize: 30,
              minSize: 30,
            },
            {
              accessorKey: 'abTarget',
              header: 'Antibody Target',
              size: 30,
              maxSize: 30,
              minSize: 30,
            },
      ],
      []
    )

    const CellTypeCart = useMemo<MRT_ColumnDef<CellTypeCart>[]>(
      () => [
          {
              accessorKey: "idCL",
              header: "Cell Ontology ID",
              size: 30,
              maxSize: 30,
              minSize: 30,
          },
          {
              accessorKey: "label",
              header: "Cell Type",
              size: 30,
              maxSize: 30,
              minSize: 30,
          },
      ],
      []
    )
    
    // Retrieve selected radio value. Afterwards, clear any existing columns
    const handleRadioSelection = (e: any) => {
      // setSelectedValue(e.target.value)
      setRadioValue(e.target.value)

      // Clear any existing data in ModelTable
      setModelColumns([])
      setModelResults([]) // Cell Types
      setModelAbKeys([]) // Antibodies
      setModelResultsDict({})// Antibodies

      // Clear any data related to plotting
      setModelPlotResults([])
      setModelPlotResultsDict(new Map())
     }

    const handleTextInputChange = (e: any) => {
      setTextInput(e.target.value);
      setSearchQuery(e.target.value);
    }

    // API call for antibody search
    const whichAntibodies = async (antibodyPayload: any) => {
        const whichAntibodiesRes = await axios.post(`${URL}/api/whichantibodies`, antibodyPayload)
        const whichAntibodiesData = whichAntibodiesRes.data
        // Renaming the "idExperiment_used" to "idExperiment" since it messes up the columns
        // Rename the key "idExperiment_used" to "idExperiment"
        whichAntibodiesData["idExperiment"] = whichAntibodiesData["idExperiment_used"];
        delete whichAntibodiesData["idExperiment_used"];
        return whichAntibodiesData
    }

    // API call for cell type search
    const whichCellTypes = async (cellTypesPayload: any) => {
        const whichCellTypesRes = await axios.post(`${URL}/api/whichcelltypes`, cellTypesPayload)
        const whichCellTypesData = whichCellTypesRes.data
        whichCellTypesData["idExperiment"] = whichCellTypesData["idExperiment_used"];
        delete whichCellTypesData["idExperiment_used"];
        return whichCellTypesData
    }

    const handleSearchButton = async (e: any) => {
        e.preventDefault()

        // Prevent action if pressing "Enter" despite button being disabled by adding extra check
        if (runDisabled) {
          return;
        }

        // Check for valid query input
        if (textInput.trim().length === 0) {
            setErrorField(true)
        } 
        
        // Handle antibody queries
        if (errorField === false && radioValue === "Antibodies") {
            const antibodyPayload = {
                "search_query": `${textInput.trim()}`
            }
            try {
                // On API call, disable the Search button
                setRunDisabled(true)

                if (textInput.trim().length !== 0) {
                    // Begin loading animation
                    setSearchIsLoading(true)
                }
                
                const antibodiesData = await whichAntibodies(antibodyPayload)
                
                if (typeof antibodiesData !== "string") {
                    // If results are not an error message, convert it into an array object
                    const antibodiesArray = convertToArrayObject(antibodiesData)
                    const antibodiesArrayID = includeID(antibodiesArray, "id")
                    setSearchResults(antibodiesArrayID)
                } else {
                    // If error string, just return an empty array
                    setSearchResults([])
                }
            } catch (error) {
                console.log(error)
            } finally {
              // Re-enable run button
              setRunDisabled(false)

              // Stop loading animation
              setSearchIsLoading(false)
            }
            
          
          // Handle cell type queries
        } else if (errorField === false && radioValue === "Cell Types") {
            const cellTypePayload = {
                "search_query": `${textInput.trim()}`
            }

            try {
                // On API call, disable the Search button
                setRunDisabled(true)

                if (textInput.trim().length !== 0) {
                    setSearchIsLoading(true)
                }

                const cellTypesData = await whichCellTypes(cellTypePayload)
         
                if (typeof cellTypesData !== "string") {
                    const cellTypesArray = convertToArrayObject(cellTypesData)
                    const cellTypesArrayID = includeID(cellTypesArray, "id")
                    setSearchResults(cellTypesArrayID)
                } else {
                    setSearchResults([])
                }
            } catch (error) {
                console.log(error)
            } finally {
              // Re-enable run button
              setRunDisabled(false)

              // Stop loading animation
              setSearchIsLoading(false)
            }
        }
        return;
    }

    // On initial load, display antibodies table and cart
    useEffect(() => {
        setSearchColumns(AntibodyColumns)
        setCartColumns(AntibodyCart)
        setRunDisabled(false) // Ensure button is disabled initially
    }, [])

    // Clear any text and search table results upon load/reload
    useEffect(() => {
        setTextInput("")
        setSearchQuery("")
        setSearchResults([])
        setSearchIsLoading(false)
        if (radioValue === "Antibodies") {
            setSearchColumns(AntibodyColumns)
            setCartColumns(AntibodyCart)
            
          } else if (radioValue === "Cell Types") {
            setSearchColumns(CellTypeColumns)
            setCartColumns(CellTypeCart)
          }
     }, [radioValue]) //dependency array, useEffect will execute whatever when this array changes

    useEffect(() => {
      setTextInput("")
    }, [radioValue])

    // Disable error message when not needed (initial load, reloads)
    useEffect(() => {
        setErrorField(false)
    }, [textInput])
 
     return (
        <div style={{display: 'flex', backgroundColor: 'white'}}>
          <div>
          <Card>
            <CardContent>
            <Typography gutterBottom sx={{ color: 'black', fontSize: 18 }}>
            Search For... </Typography>
            <FormControl>
                <RadioGroup
                  aria-labelledby="demo-radio-buttons-group-label"
                  defaultValue="Antibodies"
                  name="radio-buttons-group"
                  onChange={handleRadioSelection}
                >
                  <FormControlLabel value="Antibodies" control={<Radio />} label="Antibodies" />
                  <FormControlLabel value="Cell Types" control={<Radio />} label="Cell Types" />
                </RadioGroup>
            </FormControl>
            </CardContent>
          </Card>
          </div>

          <div style={{marginLeft: "1em", marginTop: "3em"}}>
          
          <TextField 
            error={errorField}
            helperText= {errorField == true ? "Search query cannot be blank" : ""}
            id="outlined-basic" 
            placeholder = {radioValue == "Antibodies" ? "CD4, Human, IgM..." : "Monocyte, B cell, Myeloid..."}
            variant="outlined"
            color="primary"
            onChange={handleTextInputChange}
            value={textInput} // we add this in, so that when selectedValue changes, the useEffect can see it, and reset it
            onKeyDown={(e) => e.key == "Enter" && textInput.length !== 0 ? handleSearchButton(e):null}
          />
          <Button 
            sx={{ mt: 1, mr: 1, ml: 1}} 
            type="submit" 
            variant="contained" 
            disabled={!textInput || runDisabled} 
            onClick={handleSearchButton}>
                Search
          </Button>
          </div>

        </div>
     )
}

export default RadioSearchField