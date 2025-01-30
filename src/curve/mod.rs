pub const INFINITY: f64 = f64::INFINITY; // +Inf_f64
use std::io;
use std::io::BufRead;
use std::collections::binary_heap::Iter;
use std::error::Error;
use csv::ReaderBuilder;
use csv::WriterBuilder;

pub struct Curve {
    pub m_oob_right: f64,
    pub m_oob_left: f64,
    pub m_units_x: String,
    pub m_units_y: String,
    pub m_curve_x: Vec<f64>,
    pub m_curve_y: Vec<f64>,
}

pub enum curve_axis {
    x_axis,
    y_axis,
}

pub trait CurveOperations { // un solo trait con todas las funciones o mejor varios traits (uno para cada función).
    fn read_vec(&self) -> Vec<f64>;
    fn is_oob(&self, x: f64) -> bool;
    fn set_units(&mut self, axis: curve_axis, units: String);
    fn integral(&self) -> f64;
    fn dist_mean(&self) -> f64;
    fn integrate(&mut self, k: f64) -> ();
    fn flip(&mut self) -> ();
    fn flip_2(&mut self);
    fn extend_left(&mut self) -> ();
    fn extend_right(&mut self) -> ();
    fn scale_axis_factor(&mut self, axis: curve_axis, factor: f64) -> ();
    fn scale_axis_curve(&mut self, axis: curve_axis, other: &Curve) -> ();
    fn invert_axis(&mut self, axis: curve_axis, factor: f64) -> ();
    fn x_points(&self) -> Vec<f64>;
    fn set(&mut self, x: f64, y: f64) -> ();
    fn get_point(&self, x: f64) -> f64;
    fn get_diff(&self, x: f64) -> f64;
    fn multiply_by(&mut self, other: &mut Curve) -> ();
    fn add_curve(&mut self, other: &mut Curve) -> ();
    fn add_value(&mut self, val: f64) -> ();
    fn assign(&mut self, other: &Curve) -> ();
    fn from_existing(&mut self, other: &Curve, y_units: f64) -> ();
    fn clear(&mut self) -> ();
    fn debug(&self) -> ();
    fn load(&mut self, path: &str) -> Result<(Vec<f64> , Vec<f64>), Box<dyn Error>>;
    fn save(&self, path: &str) -> Result<(), Box<dyn Error>>;
}

impl CurveOperations for Curve {

    fn read_vec(&self) -> Vec<f64> // para poder poner un vector rápido desde la terminal, el separador son los espacios.
    where // constraints on the types used in the function
        f64: std::str::FromStr, // ensures f64 can be parsed from a string using the FromStr trait
        <f64 as std::str::FromStr>::Err: std::fmt::Debug, // any error resulting from parsing f64 implements the Debug trait so it can be safely printed
    {
        std::io::stdin() // Accesses the standard input (stdin)
            .lock() // Locks the stdin for safe access, especially in multithreaded environments
            .lines() // Returns an iterator over the lines of input
            .next() // Retrieves the first line of input as an Option<Result<String>>
            .unwrap() // unwrap the Option, ensuring there is input (it panics if there's no input)
            .unwrap() // unwrap the Result, ensuring the line was successfully read (it panics if an error occurs)
            .trim() // Removes leading and trailing whitespace from the input line
            .split_whitespace() // Splits the trimmed line into substrings separated by whitespace, returning an iterator of string slices (&str)
            .map(|s| s.parse::<f64>().unwrap()) // For each substring, attempts to parse it as an f64 using s.parse::<f64>()
            .collect::<Vec<f64>>() // Collects the parsed f64 values into a Vec<f64> and returns it
    }

    fn is_oob(&self, x: f64) -> bool {
        if self.m_curve_x.is_empty() || self.m_curve_y.is_empty() {
            return true;
        }
        if x < *self.m_curve_x.first().unwrap() {
            return true;
        }
        if self.m_curve_x[self.m_curve_x.len() - 1] < x {
            return true;
        }
        else {
            return false;
        }
    }

    fn set_units(&mut self, axis: curve_axis, units: String) -> () { // Revisar esta parte del código ¿pone las mismas unidades a los dos ejes?
        match axis {
            curve_axis::x_axis => self.m_units_x = Some(units).unwrap(),
            curve_axis::y_axis => self.m_units_y = Some(units).unwrap(),
        }
    }

    fn integral(&self) -> f64 {  
        let mut accum = 0.0;
        let mut err = 0.0;
        
        let mut it = self.m_curve_x.iter();
        let mut next = it.next();

        // Out of bounds
        if self.m_oob_right + self.m_oob_left != 0.0 {
            return f64::INFINITY * (self.m_oob_left + self.m_oob_right);
        } 

        // Curva vacía
        let it = self.m_curve_x.first();
        if it == self.m_curve_x.last() {
            return 0.0;
        } 

        // normal
        else {
            //let mut index = 0;
            //while index < self.m_curve_x.len() - 1 {
            for index in 0..self.m_curve_x.len() - 1 { 
                let dx = self.m_curve_x[index+1] - self.m_curve_x[index];
                let my = 0.5 * (self.m_curve_y[index+1] + self.m_curve_y[index]);

                let term_real = my * dx - err;
                let tmp = accum + term_real;
                let term_err = tmp-accum;

                err = term_err - term_real;
                accum = tmp;

                //index = index + 1
            }
            return accum;
        }
    }

    fn dist_mean(&self) -> f64 {
        let mut accum = 0.0;
        let mut err = 0.0;
        
        let mut it = self.m_curve_x.iter();
        let mut next = it.next();

        // Out of bounds
        if self.m_oob_right + self.m_oob_left != 0.0 {
            return 0.5 * (self.m_oob_left + self.m_oob_right);
        } 

        // Curva vacía
        let it = self.m_curve_x.first();
        if it == self.m_curve_x.last() {
            return self.m_curve_x[0];
        } 

        else {
            //let mut index = 0;
            //while index < self.m_curve_x.len() - 1 {
            for index in 1..self.m_curve_x.len() - 1 {
                let dx = self.m_curve_x[index+1] - self.m_curve_x[index];
                let x0 = 0.5 * (self.m_curve_x[index+1] + self.m_curve_x[index]);
                let my = 0.5 * (self.m_curve_y[index+1] + self.m_curve_y[index]);

                let term_real = x0 * my * dx - err;
                let tmp = accum + term_real;
                let term_err = tmp-accum;

                err = term_err - term_real;
                accum = tmp;

                //index = index + 1
            }
            return accum;           
        } 
    }

    fn integrate(&mut self, k: f64) -> () {  // A RUST PARECE QUE NO LE GUSTAN LAS MAYUSCULAS

        let mut accum = 0.0;
        let mut err = 0.0;

        self.m_oob_left = k;

        if self.m_curve_x.is_empty() || self.m_curve_y.is_empty() {  // creo que la segunda condicion sobra
            self.m_oob_right = k;
            return;
        }

        if self.m_curve_x.len() == 1 {
            self.m_curve_y[0] = k;
            self.m_oob_right = k;
            return;
        }

        else {
            // Get point at x0
            let mut x0 = self.m_curve_x[0];
            let mut dx = self.m_curve_x[1] - x0;

            // Initialize the first point
            let mut x = self.m_curve_x[0];
            let mut y = self.m_curve_y[0];

            for i in 1..self.m_curve_x.len() {
                let mut x_prev = x;
                let mut y_prev = y;

                x = self.m_curve_x[i];
                y = self.m_curve_y[i];
                dx = x - x_prev;

                // Kohan sumation to calculate y:
                let term_real = 0.5 * (y + y_prev) * dx - err;
                let tmp = accum + term_real;
                let term_err = tmp - accum;
                err = term_err - term_real;
                accum = tmp;

                self.m_curve_y[i] = accum; // update curve
            }

            self.m_curve_y[0] = k; // Set the initial value at x0
            let last_dx = self.m_curve_x[self.m_curve_x.len() - 1] - self.m_curve_x[self.m_curve_x.len() - 2]; 
            let last_x_prev = self.m_curve_x[self.m_curve_x.len() - 2];
            self.m_curve_y.push(y * (x - last_x_prev + last_dx) + self.m_curve_y[self.m_curve_y.len() - 2]);

            self.m_oob_right = accum;
        }
    }

    fn flip(&mut self) -> ()  { // NO HE TERMINADO ESTA FUNCIÓN PORQUE CREO QUE HE ENCONTRADO UNA MANERA MÁS FÁCIL DE HACER EL FLIP -> VER "flip_2"
        let tmp = self.m_units_x.clone();
        self.m_units_x = self.m_units_y.clone(); // creo que no me deja mover el ownership porque la structure es mutable, lo que entra en contra con las funciones de antes -> hago .clone(), espero que sea eficiente
        self.m_units_y = tmp; 

        // ...

    }

    // A caso se podría hacer así???????????:
    fn flip_2(&mut self) {
        // Unit flip
        std::mem::swap(&mut self.m_units_x, &mut self.m_units_y);
        // x & y coordinates flip
        std::mem::swap(&mut self.m_curve_x, &mut self.m_curve_y);
    }

    fn extend_left(&mut self) -> () {
        if self.m_curve_x.is_empty() || self.m_curve_y.is_empty() {
            return;
        }
        self.m_oob_left = *self.m_curve_y.first().unwrap(); //* for deferencing the borrow (para obtener el valor real, no una referencia al valor); .unwrap() Devuelve el valor contenido en el Option, ya que first() devuelve un Option<&f64>
    }

    fn extend_right(&mut self) -> () {
        if self.m_curve_x.is_empty() || self.m_curve_y.is_empty() {
            return;
        }
        self.m_oob_right = *self.m_curve_y.last().unwrap(); // REVISAR ESTA LINEA TANTO EN EXTEND_LEFT COMO EN EXTEND_RIGHT, NO SE MUY BIEN COMO FUNCIONA, TAMPOCO LO SE EN C++
    }

    fn scale_axis_factor(&mut self, axis: curve_axis, factor: f64) -> () { // voy a probar a hacer este código con "match" ya que axis es un enum.
        match axis {
            curve_axis::x_axis => { // escalado del eje x por el factor
                for mut x in &mut self.m_curve_x {
                    *x *= factor; //* -> dereferencia: para modificar el valor al que x está apuntando (de m_curve_x), permitiendo el acceso al valor y realizar operaciones
                }
            }
            curve_axis::y_axis => { // escalado del eje y por el factor
                for mut y in &mut self.m_curve_y {
                    *y *= factor;
                }
            }
        }
        self.m_oob_left *= factor;
        self.m_oob_right *= factor;
    }

    fn scale_axis_curve(&mut self, axis: curve_axis, other: &Curve) -> () { // lo hago así para ver que tal los enums y porque yo no estoy utilizando un mapa sino dos vectores
        match axis {
            curve_axis::x_axis => {
                for i in 0..self.m_curve_x.len() {
                    self.m_curve_x[i] = other.m_curve_y[i] // REVISAR ESTO!! ME CHIRRIA!!!!!!!!!!!!!!!!!!!!
                }

            }
            curve_axis::y_axis => {
                for i in 0..self.m_curve_x.len() {
                    self.m_curve_y[i] = other.m_curve_y[i]
                }
            }
        }
        self.m_oob_left = other.m_oob_left;
        self.m_oob_right = other.m_oob_right;
    }

    fn invert_axis(&mut self, axis: curve_axis, factor: f64) -> () {
        match axis {
            curve_axis::x_axis => { // escalado del eje x por el factor
                for mut x in &mut self.m_curve_x {
                    *x = factor / *x;
                }
            }
            curve_axis::y_axis => { // escalado del eje y por el factor
                for mut y in &mut self.m_curve_y {
                    *y *= factor / *y;
                }
            }
        }
        self.m_oob_left = 1.0 / self.m_oob_left; // NO SERIA FACTOR/M_OOB_LEFT????
        self.m_oob_right = 1.0 / self.m_oob_right;
    }

    fn x_points(&self) -> Vec<f64> {
        return self.m_curve_x.clone();
    }

    fn set(&mut self, x: f64, y: f64) -> () { // ¿¿INSERTAR UN VALOR O MODIFICAR UN VALOR??
        let index = self.m_curve_x.iter().position(|&xi| xi > x).unwrap(); // para el indice donde se meterá el punto

        self.m_curve_x.insert(index, x);
        self.m_curve_y.insert(index, y);
    }

    fn get_point(&self, x: f64) -> f64 { // lo llamo get_point para que no entre en colision con otras funciones de rust
        let next = self.m_curve_x.iter().position(|&xi| xi >= x).unwrap();
        if self.m_curve_x.contains(&x) { // if there is x -> give the value
            return self.m_curve_y[next]
        }
        else { // if its not, interpolate
            let x0 = self.m_curve_x[next-1];
            let y0 = self.m_curve_x[next-1];
            let x1 = self.m_curve_x[next];
            let y1 = self.m_curve_x[next];

            let y = y0 + (x - x0) * (y1 - y0) / (x1 - x0);
            return  y;
        }
    }

    fn get_diff(&self, x: f64) -> f64 {
        let next = self.m_curve_x.iter().position(|&xi| xi >= x).unwrap();
        if self.m_curve_x[next] == *self.m_curve_x.last().unwrap() {
            return 0.0;
        }
        if self.m_curve_x[next] == *self.m_curve_x.first().unwrap() {
            return 0.0;
        }
        else {
            let x0 = self.m_curve_x[next-1];
            let y0 = self.m_curve_y[next-1];
            let x1 = self.m_curve_x[next];
            let y1 = self.m_curve_y[next];

            // two cases
            if self.m_curve_x[next] != x {
                //Middle of segment
                let diff = (y1 - y0) / (x1 - x0);
                return diff;
            }
            else {
                //Edge of a segment
                if self.m_curve_x[next+1] == *self.m_curve_x.last().unwrap() {
                    return 0.0;
                }
                else {
                    let x2 = self.m_curve_x[next+1];
                    let y2 = self.m_curve_y[next+1];

                    let diff = (y2 - y0) / (x2 - x0);
                    return diff;
                }
            }
        }
    }

    fn multiply_by(&mut self, other: &mut Curve) -> () {  // use the "get" function to obtain the points
        //combine the set of x coordinates:
        self.m_curve_x.append(&mut other.m_curve_x);
        self.m_curve_x.sort_by(|a, b| a.partial_cmp(b).unwrap()); // Como es un vector f64 no se puede usar el método tradicional .sort(), por eso hago esto que he visto en internet
        self.m_curve_x.dedup();
        // Calculo de la multiplicación
        self.m_curve_y.append(&mut other.m_curve_y); // para tamaño correcto
        for i in 0..self.m_curve_x.len() {
            self.m_curve_y[i] = self.get_point(self.m_curve_x[i]) * other.get_point(other.m_curve_x[i]);
        }
        self.m_oob_left *= other.m_oob_left;
        self.m_oob_right *= other.m_oob_right;
    }

    fn add_curve(&mut self, other: &mut Curve) -> () {     // igual que antes pero sumando 
        //combine the set of x coordinates:
        self.m_curve_x.append(&mut other.m_curve_x);
        self.m_curve_x.sort_by(|a, b| a.partial_cmp(b).unwrap()); // Como es un vector f64 no se puede usar el método tradicional .sort(), por eso hago esto que he visto en internet
        self.m_curve_x.dedup();
        // Calculo de la multiplicación
        self.m_curve_y.append(&mut other.m_curve_y); // para tamaño correcto
        for i in 0..self.m_curve_x.len() {
            self.m_curve_y[i] = self.get_point(self.m_curve_x[i]) + other.get_point(other.m_curve_x[i]);
        }
        self.m_oob_left += other.m_oob_left;
        self.m_oob_right += other.m_oob_right;
    }

    fn add_value(&mut self, val: f64) -> () {
        for i in 0..self.m_curve_x.len() {
            self.m_curve_y[i] += val;
        }
        self.m_oob_left += val;
        self.m_oob_right += val;
    }

    fn assign(&mut self, other: &Curve) -> () {
        let own_first = self.m_curve_x.first();
        let crv_first = other.m_curve_x.first();
        let own_last = self.m_curve_x.last();
        let crv_last = other.m_curve_x.last();

        // nothing to add
        if crv_first == other.m_curve_x.last() {
            return;
        }

        // no curve
        if own_first == self.m_curve_x.last() {
            return;
        }

        // Assign the middle part (must go first)
        for i in 0..self.m_curve_x.len() {
            if !other.is_oob(self.m_curve_x[i]) {
                self.m_curve_y[i] = other.get_point(self.m_curve_x[i])
            }
        }

        // if curve is longer to the left
        if crv_first < own_first {
            self.m_oob_left = other.m_oob_left;
        }

        // if curve is longer to the right 
        if own_last < crv_last {
            self.m_oob_right = other.m_oob_right;
        }

        // Assign the whole curve
        for i in 0..other.m_curve_x.len() {
            self.set(other.m_curve_x[i],other.m_curve_y[i])
        }
    }

    fn from_existing(&mut self, other: &Curve, y_units: f64) -> () {
        self.m_curve_x = other.m_curve_x.clone();
        self.m_curve_y = other.m_curve_y.clone();

        if y_units != 1.0 {
            for y in self.m_curve_y.iter_mut() {
                *y *= y_units;
            }
        }
    }

    fn clear(&mut self) -> () {
        self.m_curve_x.clear();
        self.m_curve_y.clear();
        self.m_oob_left = 0.0; // no se pueden juntar en una sola linea :(
        self.m_oob_right = 0.0;
    }

    fn debug(&self) -> () {
        for i in 0..self.m_curve_x.len() {
            print!("{} = {}, ", self.m_curve_x[i], self.m_curve_y[i])
        }
        println!();
    }

    // He dejado estas funciones para el final porque me da la sensación de que son bastante diferentes.
    // Las hago un poco a mi bola en vez de intentar replicar las de curve.cpp

    fn load(&mut self, path: &str) -> Result<(Vec<f64> , Vec<f64>), Box<dyn Error>> { // De momento solo cargo las coordenadas de la curva, sencillito

        let mut rdr = ReaderBuilder::new().has_headers(false).from_path(path)?;  //NO SE QUE SIGNIFICA "?", creo que es para propagación de errores?
        for result in rdr.records() {
            let record = result?;

            let x: f64 = record.get(0).unwrap().parse()?;
            let y: f64 = record.get(1).unwrap().parse()?;

            self.m_curve_x.push(x);
            self.m_curve_y.push(y);
        }
        Ok((self.m_curve_x.clone(),self.m_curve_y.clone())) //para el return
    }

    fn save(&self, path: &str) -> Result<(), Box<dyn Error>> {
        
        let mut wtr = WriterBuilder::new().has_headers(false).from_path(path).unwrap();
        //wtr.write_record(&["X", "Y"])?; // headers, para tener de referencia por si acaso
        for (x,y) in self.m_curve_x.iter().zip(self.m_curve_y.iter()) {
            wtr.write_record(&[x.to_string(),y.to_string()]).unwrap();
        }
        wtr.flush().unwrap(); //no se por que se hace esto pero en los ejemplos siempre lo hacen
        Ok(())
    }
}