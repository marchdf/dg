// Commentaires

lc = 0.3;   // longueur caractéristique

// Définition des Points

Point(1) = {-1,0,0,lc};
Point(2) = {1,0,0,lc};

// Définition des Lignes

Line(1) = {1,2};

// Définition de la Surface

Physical Line(1) = {1,2,3,4};

Transfinite Line{1}=NUM;

