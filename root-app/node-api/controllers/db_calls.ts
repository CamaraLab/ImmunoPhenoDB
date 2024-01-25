import axios from "axios"
import { Request, Response} from "express"
import dotenv from "dotenv";

dotenv.config({ path: "../.env"})
const FLASK_PORT = 5000

// Deploying onto the cloud
const cloud_docker_URL = `http://flask:${FLASK_PORT}`

const system_URL = cloud_docker_URL

console.log("FLASK URL:", system_URL)

declare module "express-serve-static-core" {
    interface Request {
        [key: string]: Array<string>;
    }
}

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