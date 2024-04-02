import axios from "axios"
import { createClient } from 'redis';
import { Request, Response} from "express"
import dotenv from "dotenv";

dotenv.config({ path: "../.env"})
const FLASK_PORT = 5000

// Deploying onto the cloud
const system_URL = `http://flask:${FLASK_PORT}`

console.log("FLASK URL:", system_URL)

declare module "express-serve-static-core" {
    interface Request {
        [key: string]: Array<string>;
    }
}

const client: ReturnType<typeof createClient> = createClient({
    url: 'redis://redis:6379' 
});
client.on('error', err => console.log('Redis Client Error', err));
console.log("trying to connect to redis...")
client.connect();
console.log("redis connected")

export const idCLs = async (req: any, res: Response) => {
    try {
        const idCLs_list = await axios.get(`${system_URL}/idcls`)
        console.log("unique idCLs:", idCLs_list)
        res.status(200).send(idCLs_list.data)
    } 
    catch (error) {
        if (error) {
            console.log(error)
            res.status(400).send("Bad request. Error with idCLs server request.")
        }
    }
}

export const tissues = async(req: any, res: Response) => {
    try{
        const flask_tissues = await axios.get(`${system_URL}/tissues`)
        console.log("tissues:", flask_tissues)
        res.status(200).send(flask_tissues.data)
    }
    catch (error) {
        if (error) {
            console.log(error)
            res.status(400).send("Bad request. Errors with tissues server request.")
        }
    }
}

export const findAbs = async (req: Request, res: Response) => {
    try {
        console.log("findAbs recieved from client:", req.body)
        const flask_abs = await axios.post(`${system_URL}/findabs`, req.body)
        console.log("findAbs final result:", flask_abs)
        res.status(200).send(flask_abs.data)
    }
    catch (error) {
        if (error) {
            console.log(error)
            res.status(400).send("Bad request. Error with findAbs computation.")
        }
    }
}

export const plotAbs = async (req: any, res: any) => {
    try {
        console.log("plotAbs recieved from client", req.body)
        const flask_plotabs = await axios.post(`${system_URL}/plotabs`, req.body)
        // console.log("plotAbs final result:", flask_plotabs)
        res.status(200).send(flask_plotabs.data)
    }
    catch (error) {
        if (error) {
            console.log(error)
            res.status(400).send("Bad request. Error with plotAbs computation.")
        }
    }
}

export const findCellTypes = async (req: any, res: any) => {
    try {
        console.log("findCellTypes recieved from client:", req.body)
        const flask_celltypes = await axios.post(`${system_URL}/findcelltypes`, req.body)
        console.log("findCellTypes final result:", flask_celltypes)
        res.status(200).send(flask_celltypes.data)
    }
    catch (error) {
        if (error) {
            console.log(error)
            res.status(400).send("Bad request. Error with findCellTypes computation.")
        }
    }
}

export const plotCellTypes = async (req: any, res: any) => {
    try {
        console.log("plotCellTypes recieved from client", req.body)
        const flask_plotcelltypes = await axios.post(`${system_URL}/plotcelltypes`, req.body)
        console.log("plotCellTypes final result:", flask_plotcelltypes)
        res.status(200).send(flask_plotcelltypes.data)
    }
    catch (error) {
        if (error) {
            console.log(error)
            res.status(400).send("Bad request. Error with plotCellTypes computation.")
        }
    }
}

export const findExperiments = async (req: any, res: any) => {
    try {
        console.log("getExperiments recieved from client", req.body)
        const flask_findexperiments = await axios.post(`${system_URL}/findexperiments`, req.body)
        console.log("getExperiments final result:", flask_findexperiments)
        res.status(200).send(flask_findexperiments.data)
    }
    catch (error) {
        if (error) {
            console.log(error)
            res.status(400).send("Bad request. Error with findExperiments computation.")
        }
    }
}

export const whichAntibodies = async (req: any, res: any) => {
    try {
        console.log("whichAntibodies recieved from client", req.body)
        const flask_whichantibodies = await axios.post(`${system_URL}/whichantibodies`, req.body)
        console.log("which antibodies:", flask_whichantibodies)
        res.status(200).send(flask_whichantibodies.data)
    } 
    catch (error) {
        if (error) {
            console.log(error)
            res.status(400).send("Bad request. Error with whichAntibodies computation.")
        }
    }
}

export const whichCellTypes = async (req: any, res: any) => {
    try {
        console.log("whichCellTypes recieved from client", req.body)
        const flask_whichcelltypes = await axios.post(`${system_URL}/whichcelltypes`, req.body)
        console.log("which celltypes:", flask_whichcelltypes)
        res.status(200).send(flask_whichcelltypes.data)
    } 
    catch (error) {
        if (error) {
            console.log(error)
            res.status(400).send("Bad request. Error with whichCellTypes computation.")
        }
    }
}

export const whichExperiments = async (req: any, res: any) => {
    try{
        console.log("whichExperiments recieved from client", req.body)
        const flask_whichexperiments = await axios.post(`${system_URL}/whichexperiments`, req.body)
        console.log("which experiments:", flask_whichexperiments)
        res.status(200).send(flask_whichexperiments.data)
    }
    catch (error) {
        if (error) {
            console.log(error)
            res.status(400).send("Bad request. Error with whichExperiments computation.")
        }
    }
}

export const stveaReference = async (req: Request, res: Response) => {
    try {
        console.log("stveaReference recieved from client:", req.body)
        const flask_stvea_Reference = await axios.post(`${system_URL}/stveareference`, req.body)
        // console.log("stveaReference final result:", flask_stvea_Reference)
        res.status(200).send(flask_stvea_Reference.data)
    }
    catch (error) {
        if (error) {
            console.log(error)
            res.status(400).send("Bad request. Error with stveaReference computation.")
        }
    }
}

export const convertCellType = async (req: Request, res: Response) => {
    console.log("we entered convertCellType")
    async function fetchValuesFromAPI(keys: any) {
        const apiResults = [];
      
        for (const key of keys) {
            console.log("we got this missing key:", key)
          try {
            const idCL_params = {
                'q': key,
                'exact': 'true',
                'ontology': 'cl',
                'fieldList': 'label',
                'rows': 1,
                'start': 0
            }
    
            // Replace the following URL with your actual API endpoint
            const response = await axios.get('https://www.ebi.ac.uk/ols4/api/search', {params: idCL_params});
            // console.log("server response:", response)
            const value = response.data; // Adjust this based on your API response structure
            const cellType = value['response']['docs'][0]['label']
            console.log("calling external API..... value:", value)
            console.log("Extracting from value:", cellType)
            apiResults.push({ key, cellType });
    
    
          } catch (error: any) {
            console.error(`Error fetching value for key ${key}:`, error.message);
          }
        }
        return apiResults;
      }

      try {
        const keys = req.body["idCL"]; // Assuming the request body is an array of keys ["A", "B", "C"]
        // console.log("body:",  req.body);
        console.log("fetching redis for these keys:", keys)
        const results: any = {};
    
        // Step 2: Use mGet to fetch values from Redis
        const mGetRes = await client.mGet(keys)
        console.log("redis results mGetRes:", mGetRes)
  
        // Step 3: Process values from Redis
        mGetRes.forEach((value, index) => {
          if (value != null) {
            //   console.log("entered here! valid redis result")
              // Key exists in Redis, add to results
              results[keys[index]] = value;
          }
        })
  
        // Step 4: Fetch values for keys not found in Redis
        const missingKeys = keys.filter((key: string, index: number) => mGetRes[index] === null);
        console.log("missingKeys:", missingKeys);
        const apiResults = await fetchValuesFromAPI(missingKeys);
        console.log("apiResults for missing values:", apiResults)
  
        // Step 5: Store new values in Redis and add to results
        apiResults.forEach(({ key, cellType }) => {
            client.set(key, cellType);
            results[key] = cellType;
        });
  
        // Print final results dict
        console.log("results dict:", results);
  
        // Step 6: send final dict back to client
        res.status(200).json({ success: true, results});
       
      } 
      catch (error) {
        console.error('Error processing keys:', error);
        res.status(400).json("Bad request. Error with convertCellType computation.");
      }

}