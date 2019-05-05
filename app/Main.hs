module Main where

import VectorSpace

main :: IO ()
main =
    putStrLn $ show $
    orthProject
        (IntegralInnerProduct 0 1)
        (RealFunction exp)
        ([Polynomial [1, 0, 0, 0], Polynomial [0, 1, 0, 0], Polynomial [0, 0, 1, 0], Polynomial [0, 0, 0, 1]] :: [Polynomial 4 Double])
