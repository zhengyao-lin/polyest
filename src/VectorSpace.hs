{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE QuantifiedConstraints #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GADTs #-}

{- basic definitions for fields and vector spaces -}

module VectorSpace where

import GHC.TypeLits

import Data.List
import Data.Proxy

import Numeric.Integration.TanhSinh

import Debug.Trace

infinity = read "Infinity"

class Num f => Field f where

class Fractional f => InnerProductField f where
    sqrt :: f -> f

-- basic fields: real(double) and rational
instance Field Double
instance Field Rational

instance InnerProductField Double where
    sqrt = (** 0.5)

class Field f => VectorSpace s f where
    infixl 6 +.
    infixl 6 -.
    infixl 7 *.

    vs :: s f
    vs = error "cannot be evaluated"

    zero :: s f
    
    dim :: Integral t => s f -> Maybe t

    (+.) :: s f -> s f -> s f
    (*.) :: f -> s f -> s f

    (-.) :: s f -> s f -> s f
    v1 -. v2 = v1 +. (-1) *. v2

class (VectorSpace s f, VectorSpace s' f) => Subspace s' s f where
    super :: s' f -> s f

instance VectorSpace s f => Subspace s s f where
    super = id

class (Field f, VectorSpace s f) => FunctionSpace s f where
    eval :: s f -> f -> f

class (InnerProductField f, VectorSpace s f) => InnerProductSpace i s f where
    inner :: i -> s f -> s f -> f

    norm :: i -> s f -> f
    norm i v = VectorSpace.sqrt (inner i v v)

class FunctionSpace s Double => Integrable s where
    integrate :: s Double -> Double -> Double -> Double
    integrate f a b = result $ absolute infinity $ simpson (eval f) a b -- quadRes (quadRomberg defQuad (a, b) (eval f))

-- basic vectors by cartesian product
data Vector (n :: Nat) f = Vector [f] deriving (Eq, Show)

instance (Field f, KnownNat n) => VectorSpace (Vector n) f where
    zero = Vector (replicate d 0)
        where Just d = dim (vs :: Vector n f)

    (Vector v1) +. (Vector v2) = Vector (zipWith (+) v1 v2)
    c *. (Vector v) = Vector (map (c *) v)

    dim _ = Just (fromIntegral (natVal (Proxy :: Proxy n)))

-- polynomial space
data Polynomial (n :: Nat) c = Polynomial [c] deriving (Eq)

instance (Field f, KnownNat n) => VectorSpace (Polynomial n) f where
    zero = Polynomial (replicate d 0)
        where Just d = dim (vs :: Polynomial n f)

    (Polynomial p1) +. (Polynomial p2) =
        Polynomial (zipWith (+) p1 p2)
        where
            l1 = length p1
            l2 = length p2
            (p1', p2') =
                if l1 < l2 then
                    (p1 ++ replicate (l2 - l1) 0, p2)
                else
                    (p1, p2 ++ replicate (l1 - l2) 0)

    c *. (Polynomial p) = Polynomial (map (c *) p)

    dim _ = Just (fromIntegral (natVal (Proxy :: Proxy n)))

instance (Field f, KnownNat n) => FunctionSpace (Polynomial n) f where
    eval (Polynomial p) x = sum (zipWith (\c p -> c * x ^ p) p [0..])

instance KnownNat n => Integrable (Polynomial n)

instance KnownNat n => Subspace (Polynomial n) RealFunction Double where
    super f = RealFunction (eval f)

instance KnownNat n => Show (Polynomial n Double) where
    show (Polynomial p) = intercalate "+" (map (\(c, p) -> show c ++ "x^" ++ show p) (zip p [0..]))

-- function space restricted to some integer interval
data RealFunction f = RealFunction (f -> f)

toRealFunction :: FunctionSpace s Double => s Double -> RealFunction Double
toRealFunction f = RealFunction (eval f)

-- lowerRange :: forall a b f t. (KnownNat a, Num t) => RealFunction a b f -> t
-- lowerRange _ = fromIntegral (natVal (Proxy :: Proxy a))

-- upperRange :: forall a b f t. (KnownNat b, Num t) => RealFunction a b f -> t
-- upperRange _ = fromIntegral (natVal (Proxy :: Proxy b))

-- inRange :: (KnownNat a, KnownNat b) => RealFunction a b f -> Double -> Bool
-- inRange v x = lowerRange v <= x && x <= upperRange v

instance VectorSpace RealFunction Double where
    zero = RealFunction (const 0)

    (RealFunction f1) +. (RealFunction f2) = RealFunction (\x -> f1 x + f2 x)

    c *. (RealFunction f) = RealFunction (\x -> c * f x)

    dim _ = error "infinite dim"

instance FunctionSpace RealFunction Double where
    eval (RealFunction f) x = f x

instance Integrable RealFunction

-- define various inner products

data StandardInnerProduct = StandardInnerProduct deriving (Eq, Show)
data IntegralInnerProduct = IntegralInnerProduct Double Double deriving (Eq, Show)

instance KnownNat n => InnerProductSpace StandardInnerProduct (Vector n) Double where
    inner _ (Vector v1) (Vector v2) = sum (zipWith (*) v1 v2)

instance KnownNat n => InnerProductSpace StandardInnerProduct (Polynomial n) Double where
    inner _ (Polynomial p1) (Polynomial p2) = sum (zipWith (*) p1 p2)

instance Integrable f => InnerProductSpace IntegralInnerProduct f Double where
    inner (IntegralInnerProduct a b) f1 f2 =
        integrate (RealFunction (\x -> eval f1 x * eval f2 x)) a b

-- common utilities
gramSchmidt :: InnerProductSpace i s f => i -> [s f] -> [s f]
gramSchmidt i basis =
    foldl extend [] basis
    where
        extend basis' v =
            let e = foldl (\d v' -> d -. inner i v v' *. v') v basis'
                e' = (1 / norm i e) *. e
            in basis' ++ [e']

orthProject :: (InnerProductSpace i s f, InnerProductSpace i s' f, Subspace s' s f)
               => i -> s f -> [s' f] -> s' f
orthProject i f basis =
    foldl (+.) zero (map (\v -> inner i f (super v) *. v) basis')
    where basis' = gramSchmidt i basis -- get orth basis
