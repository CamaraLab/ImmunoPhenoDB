import express, { Express, Request, Response } from 'express'
import dotenv from 'dotenv';
import router from './routes/db_routes.js'

dotenv.config({ path: '../.env'})

const app: Express = express()
app.use(express.json())

// const NODE_PORT = process.env.NODE_PORT
const NODE_PORT = 3000
app.listen(NODE_PORT, () => console.log(`Server is up and running at port ${NODE_PORT}`))

app.use('/', router)