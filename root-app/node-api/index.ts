import express, { Express, Request, Response } from 'express'
import dotenv from 'dotenv';
import router from './routes/db_routes.js'

dotenv.config({ path: '../.env'})

import cors from 'cors';

const app: Express = express()
app.use(express.json())

app.use(cors())

// const NODE_PORT = process.env.NODE_PORT
const NODE_PORT = 8080
app.listen(NODE_PORT, () => console.log(`Node server is up and running at port ${NODE_PORT}`))

app.use('/api', router)