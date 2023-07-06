import axios from "axios"
import { Request, Response} from "express"
import dotenv from "dotenv";

dotenv.config({ path: "../.env"})
// const FLASK_PORT = process.env.FLASK_PORT
const FLASK_PORT = 5000
const flask_URL = `http://flask:${FLASK_PORT}`
console.log("FLASK URL:", flask_URL)

declare module "express-serve-static-core" {
    interface Request {
        [key: string]: Array<string>;
    }
}

export const idCLs = async (req: any, res: Response) => {
    try {
        const idCLs_list = await axios.get(`${flask_URL}/idCLs`)
        console.log("unique idCLs:", idCLs_list)
        res.status(200).send(idCLs_list.data)
    } 
    catch (error) {
        if (error) {
            console.log(error)
            res.send()
        }
    }
}

export const runLMM = async (req: Request, res: Response) => {
    try {
        console.log("recieved from client:", req.body)
        const flask_json = await axios.post(`${flask_URL}/abLMM`, req.body)
        console.log("final result:", flask_json)
        res.status(200).send(flask_json.data)
    }
    catch (error) {
        if (error) {
            console.log(error)
            res.send()
        }
    }
    
}